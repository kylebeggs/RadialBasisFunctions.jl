"""
    Regrid

Builds an operator for interpolating from one set of points to another.
"""
struct Regrid
    ℒ
    Regrid() = new(identity)
end

# convienience constructors
"""
    function regrid(data, eval_points, order, dim, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is an regrid from one set of points to another, `data` -> `eval_points`.
"""
function regrid(
    data::AbstractVector,
    eval_points::AbstractVector,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
    adjl=find_neighbors(data, eval_points, k),
) where {T<:Int,B<:AbstractRadialBasis}
    return RadialBasisOperator(Regrid(), data, eval_points, basis; k=k, adjl=adjl)
end

# pretty printing
print_op(op::Regrid) = "regrid"
