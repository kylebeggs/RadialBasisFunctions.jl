abstract type AbstractOperator end
abstract type ScalarValuedOperator <: AbstractOperator end
abstract type VectorValuedOperator <: AbstractOperator end

"""
    struct RadialBasisOperator

Operator of data using a radial basis with potential monomial augmentation.
"""
struct RadialBasisOperator{L,W,D,C,A,B<:AbstractRadialBasis}
    ℒ::L
    weights::W
    data::D
    eval_points::C
    adjl::A
    basis::B
    valid_cache::Base.RefValue{Bool}
    function RadialBasisOperator(
        ℒ::L, weights::W, data::D, eval_points::C, adjl::A, basis::B
    ) where {L,W,D,C,A,B<:AbstractRadialBasis}
        return new{L,W,D,C,A,B}(ℒ, weights, data, eval_points, adjl, basis, Ref(false))
    end
end

# convienience constructors
function RadialBasisOperator(
    ℒ, data::AbstractVector{D}, basis::B=PHS(3; poly_deg=2); k::T=autoselect_k(data, basis)
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    Na = length(adjl)
    Nd = length(data)
    weights = spzeros(eltype(D), Na, Nd)
    return RadialBasisOperator(ℒ, weights, data, data, adjl, basis)
end

function RadialBasisOperator(
    ℒ,
    data::AbstractVector{D},
    eval_points::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, eval_points, k)
    Na = length(adjl)
    Nd = length(data)
    weights = spzeros(eltype(D), Na, Nd)
    return RadialBasisOperator(ℒ, weights, data, eval_points, adjl, basis)
end

# extend Base methods
Base.length(op::RadialBasisOperator) = length(op.adjl)
Base.size(op::RadialBasisOperator) = size(op.weights)
function Base.size(op::RadialBasisOperator{<:VectorValuedOperator})
    return ntuple(i -> size(op.weights[i]), embeddim(op))
end
Base.getindex(op::RadialBasisOperator, i) = nonzeros(op.weights[i, :])
function Base.getindex(op::RadialBasisOperator{VectorValuedOperator}, i)
    return ntuple(j -> nonzeros(op.weights[j][i, :]), embeddim(op))
end

# convienience methods
embeddim(op::RadialBasisOperator) = length(first(op.data))

# caching
invalidate_cache(op::RadialBasisOperator) = op.valid_cache[] = false
validate_cache(op::RadialBasisOperator) = op.valid_cache[] = true
is_cache_valid(op::RadialBasisOperator) = op.valid_cache[]

# dispatches for evaluation
_eval_op(op::RadialBasisOperator, x::AbstractVector) = _eval_op(op.weights, x)

function _eval_op(op::RadialBasisOperator{<:VectorValuedOperator}, x::AbstractVector)
    return ntuple(i -> _eval_op(op.weights[i], x), embeddim(op))
end

function _eval_op(
    op::RadialBasisOperator{<:VectorValuedOperator,W}, x::AbstractVector
) where {W<:Tuple}
    if first(op.weights) isa Vector{<:Vector}
        return ntuple(i -> _eval_op(op.weights[i], x, op.adjl), embeddim(op))
    else
        return ntuple(i -> _eval_op(op.weights[i], x), embeddim(op))
    end
end

function _eval_op(
    op::RadialBasisOperator{L,W}, x::AbstractVector
) where {L,W<:Vector{<:Vector}}
    return _eval_op(op.weights, x, op.adjl)
end

_eval_op(w::AbstractMatrix, x::AbstractVector) = w * x

function _eval_op(w::AbstractVector{<:AbstractVector{T}}, x::AbstractVector, adjl) where {T}
    y = zeros(T, length(w))
    Threads.@threads for i in eachindex(adjl)
        @views y[i] = w[i] ⋅ x[adjl[i]]
    end
    return y
end

# evaluate
function (op::RadialBasisOperator)(x)
    !is_cache_valid(op) && update_weights!(op)
    return _eval_op(op, x)
end

function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::RadialBasisOperator, x::AbstractVecOrMat
)
    !is_cache_valid(op) && update_weights!(op)
    return mul!(y, op.weights, x)
end
function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::RadialBasisOperator, x::AbstractVecOrMat, α, β
)
    !is_cache_valid(op) && update_weights!(op)
    return mul!(y, op.weights, x, α, β)
end

function LinearAlgebra.:⋅(
    op::RadialBasisOperator{<:VectorValuedOperator}, x::AbstractVector
)
    !is_cache_valid(op) && update_weights!(op)
    return sum(op(x))
end

function LinearAlgebra.mul!(
    y::AbstractVector{<:Real},
    op::RadialBasisOperator{<:VectorValuedOperator},
    x::AbstractVector,
)
    !is_cache_valid(op) && update_weights!(op)
    for i in eachindex(op.weights)
        mul!(y[i], op.weights[i], x)
    end
end

# update weights
function update_weights!(op::RadialBasisOperator)
    op.weights .= _build_weights(op.ℒ.ℒ, op)
    validate_cache(op)
    return nothing
end

function update_weights!(op::RadialBasisOperator{<:VectorValuedOperator})
    return _update_weights!(op, op.weights)
end

function _update_weights!(op, weights::NTuple{N,AbstractMatrix}) where {N}
    for (i, ℒ) in enumerate(op.ℒ.ℒ)
        weights[i] .= _build_weights(ℒ, op)
    end
    validate_cache(op)
    return nothing
end

function _update_weights!(op, weights::NTuple{N,AbstractVector}) where {N}
    for (i, ℒ) in enumerate(op.ℒ.ℒ)
        w = _build_weights(ℒ, op)
        for j in eachindex(weights[i])
            weights[i][j] .= w[j]
        end
    end
    validate_cache(op)
    return nothing
end

# pretty printing
function Base.show(io::IO, op::RadialBasisOperator)
    println(io, "RadialBasisOperator")
    println(io, "├─Operator: " * print_op(op.ℒ))
    println(io, "├─Data type: ", typeof(first(op.data)))
    println(io, "├─Number of points: ", length(op.data))
    println(io, "├─Stencil size: ", length(first(op.adjl)))
    return println(
        io,
        "└─Basis: ",
        print_basis(op.basis),
        " with degree $(op.basis.poly_deg) polynomial augmentation",
    )
end

print_op(op) = "$op"
