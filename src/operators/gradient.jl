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
            x -> ∂(x, 1, dim)
        end
    end
    ℒ = Gradient(f)
    return RadialBasisOperator(ℒ, data, basis; k=k)
end

function RadialBasisOperator(
    ℒ::Gradient,
    data::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    N = length(data)
    weights = ntuple(_ -> spzeros(N, N), length(ℒ.ℒ))
    return RadialBasisOperator(ℒ, weights, data, adjl, basis)
end

# pretty printing
print_op(op::Gradient) = "Gradient (∇)"
