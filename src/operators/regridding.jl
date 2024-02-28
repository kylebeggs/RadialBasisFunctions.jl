"""
    Interpolator

Builds an operator for interpolating from one set of points to another.
"""
struct Interpolator{L}
    ℒ::L
    Interpolator() = new{typeof(identity)}(identity)
end

# convienience constructors
"""
    function interpolator(data, eval_points, order, dim, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the partial derivative, `Partial`. The resulting operator will only evaluate at `eval_points`.
"""
function interpolator(
    data::AbstractVector{D},
    eval_points::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    return RadialBasisOperator(Interpolator(), data, eval_points, basis; k=k)
end

# pretty printing
print_op(op::Interpolator) = "identity"
