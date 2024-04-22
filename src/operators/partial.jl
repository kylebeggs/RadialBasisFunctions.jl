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
"""
    function partial(data, order, dim, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the partial derivative, `Partial`, of `order` with respect to `dim`.
"""
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

"""
    function partial(data, eval_points, order, dim, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the partial derivative, `Partial`. The resulting operator will only evaluate at `eval_points`.
"""
function partial(
    data::AbstractVector,
    eval_points::AbstractVector,
    order::T,
    dim::T,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {T<:Int,B<:AbstractRadialBasis}
    f = let o = order, dim = dim
        x -> ∂(x, o, dim)
    end
    ℒ = Partial(f, order, dim)
    return RadialBasisOperator(ℒ, data, eval_points, basis; k=k)
end

# pretty printing
print_op(op::Partial) = "∂ⁿf/∂xᵢ (n = $(op.order), i = $(op.dim))"
