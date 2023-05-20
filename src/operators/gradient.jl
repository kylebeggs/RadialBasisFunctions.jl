"""
    Gradient <: VectorValuedOperator

Builds an operator for the gradient of a function.
"""
struct Gradient{L<:NTuple} <: VectorValuedOperator
    ℒ::L
end

# convienience constructors
function gradient(
    data::AbstractVector{D}, basis::B=PHS(3; poly_deg=2); k::T=autoselect_k(data, basis)
) where {D<:AbstractArray,B<:AbstractRadialBasis,T<:Int}
    f = ntuple(length(first(data))) do dim
        return let dim = dim
            x -> ∂(x, Val{1}(), dim)
        end
    end
    ℒ = Gradient(f)
    return RadialBasisOperator(ℒ, data, basis; k=k)
end

# evaluate
function (op::RadialBasisOperator{VectorValuedOperator})(x::AbstractVector)
    return map(w -> w * x, op.weights)
end
function LinearAlgebra.:⋅(op::RadialBasisOperator{VectorValuedOperator}, x::AbstractVector)
    return sum(op(x))
end

# TODO
function LinearAlgebra.mul!(
    y::AbstractVector{<:Real}, op::VectorValuedOperator, x::AbstractVector
)
    for i in eachindex(op.weights)
        mul!(y[i], op.weights[i], x)
    end
end

# pretty printing
print_op(op::Gradient) = "Gradient (∇)"
