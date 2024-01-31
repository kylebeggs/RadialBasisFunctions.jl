"""
    Partial <: ScalarValuedOperator

Builds an operator for a first order partial derivative.
"""
struct Partial{L<:Function,T<:Int} <: ScalarValuedOperator
    ℒ::L
    order::T
    dim::T
end

# convienience constructors
function partial(
    data::AbstractVector{D},
    order::T,
    dim::T,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    f = let o = order, dim = dim
        x -> ∂(x, o, dim)
    end
    ℒ = Partial(f, order, dim)
    return RadialBasisOperator(ℒ, data, basis; k=k)
end

function partial(
    data::AbstractVector{D},
    centers::AbstractVector{D},
    order::T,
    dim::T,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    f = let o = order, dim = dim
        x -> ∂(x, o, dim)
    end
    ℒ = Partial(f, order, dim)
    return RadialBasisOperator(ℒ, data, centers, basis; k=k)
end

# pretty printing
print_op(op::Partial) = "∂ⁿf/∂xᵢ (n = $(op.order), i = $(op.dim))"
