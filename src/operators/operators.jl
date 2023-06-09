abstract type AbstractRadialBasisOperator end
abstract type AbstractOperator end
abstract type ScalarValuedOperator <: AbstractOperator end
abstract type VectorValuedOperator <: AbstractOperator end

"""
    struct RadialBasisOperator <: AbstractRadialBasisOperator

Operator of data using a radial basis with potential monomial augmentation.
"""
struct RadialBasisOperator{L,W,D,A,B<:AbstractRadialBasis} <: AbstractRadialBasisOperator
    ℒ::L
    weights::W
    data::D
    adjl::A
    basis::B
    valid_cache::Base.RefValue{Bool}
    function RadialBasisOperator(
        ℒ::L, weights::W, data::D, adjl::A, basis::B
    ) where {L,W,D,A,B<:AbstractRadialBasis}
        return new{L,W,D,A,B}(ℒ, weights, data, adjl, basis, Ref(false))
    end
end

# convienience constructors
function RadialBasisOperator(
    ℒ, data::AbstractVector{D}, basis::B=PHS(3; poly_deg=2); k::T=autoselect_k(data, basis)
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    N = length(data)
    weights = spzeros(N, N)
    return RadialBasisOperator(ℒ, weights, data, adjl, basis)
end

# evaluate
function (op::RadialBasisOperator)(x)
    !is_cache_valid(op) && update_weights!(op)
    return op.weights * x
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

# evaluate
function (op::RadialBasisOperator{<:VectorValuedOperator})(x::AbstractVector)
    !is_cache_valid(op) && update_weights!(op)
    return map(w -> w * x, op.weights)
end
function LinearAlgebra.:⋅(
    op::RadialBasisOperator{<:VectorValuedOperator}, x::AbstractVector
)
    !is_cache_valid(op) && update_weights!(op)
    return sum(op(x))
end

# TODO
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
    op.weights .= _build_weightmx(op.ℒ, op.data, op.adjl, op.basis)
    validate_cache(op)
    return nothing
end

function update_weights!(op::RadialBasisOperator{<:VectorValuedOperator})
    for (i, ℒ) in enumerate(op.ℒ.ℒ)
        op.weights[i] .= _build_weightmx(ℒ, op.data, op.adjl, op.basis)
    end
    validate_cache(op)
    return nothing
end

Base.getindex(op::O, i) where {O<:AbstractRadialBasisOperator} = nonzeros(op.weights[i, :])
invalidate_cache(op::RadialBasisOperator) = op.valid_cache[] = false
validate_cache(op::RadialBasisOperator) = op.valid_cache[] = true
is_cache_valid(op::RadialBasisOperator) = op.valid_cache[]

(op::AbstractOperator)(x) = op.ℒ(x)

# pretty printing
function Base.show(io::IO, op::RadialBasisOperator)
    println(io, "RadialBasisOperator")
    println(io, "  └─Operator: " * print_op(op.ℒ))
    println(io, "  └─Data type: ", typeof(first(op.data)))
    println(io, "  └─Number of points: ", length(op.data))
    println(io, "  └─Dimensions: ", length(first(op.data)))
    println(io, "  └─Stencil size: ", length(first(op.adjl)))
    return println(
        io,
        "  └─Basis: ",
        print_basis(op.basis),
        " with degree $(op.basis.poly_deg) Monomial",
    )
end

print_op(op) = "$op"
