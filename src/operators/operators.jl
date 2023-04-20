"""
    struct RadialBasisOperator <: AbstractRadialBasisOperator

Operator of data using a radial basis with potential monomial augmentation.
"""
struct RadialBasisOperator{L,W<:AbstractMatrix,D,A,B<:AbstractRadialBasis} <:
       AbstractRadialBasisOperator
    ℒ::L
    weights::W
    data::D
    adjl::A
    basis::B
    valid_cache::Base.RefValue{Bool}
    function RadialBasisOperator(
        ℒ, weights::W, data::D, adjl::A, basis::B
    ) where {W<:AbstractMatrix,D,A,B<:AbstractRadialBasis}
        return new{typeof(ℒ),W,D,A,B}(ℒ, weights, data, adjl, basis, Ref(false))
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
    !op.valid_cache[] && update_weights!(op)
    return op.weights * x
end

# update weights
function update_weights!(op::RadialBasisOperator)
    op.weights .= _build_weightmx(op.ℒ, op.data, op.adjl, op.basis)
    op.valid_cache[] = true
    return nothing
end

# operations
function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::RadialBasisOperator, x::AbstractVecOrMat
)
    !op.valid_cache[] && update_weights!(op)
    return mul!(y, op.weights, x)
end
function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::RadialBasisOperator, x::AbstractVecOrMat, α, β
)
    !op.valid_cache[] && update_weights!(op)
    return mul!(y, op.weights, x, α, β)
end

Base.getindex(op::O, i) where {O<:AbstractRadialBasisOperator} = nonzeros(op.weights[i, :])

# include built-in operators
include("scalar_valued.jl")
include("vector_valued.jl")
include("operator_operations.jl")

# pretty printing
function Base.show(io::IO, op::RadialBasisOperator)
    println(io, "RadialBasisOperator")
    println(io, "  Operator: " * print_op(op.ℒ))
    println(io, "  Data type: ", typeof(first(op.data)))
    println(io, "  Number of points: ", length(op.data))
    println(io, "  Dimensions: ", length(first(op.data)))
    println(io, "  Stencil size: ", length(first(op.adjl)))
    return println(
        io, "  Basis: ", typeof(op.basis), " with degree $(op.basis.poly_deg) Monomial"
    )
end

print_op(op) = "$op"
