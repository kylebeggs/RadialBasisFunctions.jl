"""
    Directional <: ScalarValuedOperator

Operator for the directional derivative, or the inner product of the gradient and a direction vector.
"""
struct Directional{L<:NTuple,T} <: ScalarValuedOperator
    ℒ::L
    v::T
end

"""
    function directional(data, v, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the directional derivative, `Directional`.
"""
function directional(
    data::AbstractVector{D},
    v::AbstractVector,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,B<:AbstractRadialBasis,T<:Int}
    f = ntuple(length(first(data))) do dim
        return let dim = dim
            x -> ∂(x, 1, dim)
        end
    end
    ℒ = Directional(f, v)
    return RadialBasisOperator(ℒ, data, basis; k=k)
end

"""
    function directional(data, eval_points, v, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the directional derivative, `Directional`.
"""
function directional(
    data::AbstractVector,
    eval_points::AbstractVector,
    v::AbstractVector,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {B<:AbstractRadialBasis,T<:Int}
    f = ntuple(length(first(data))) do dim
        return let dim = dim
            x -> ∂(x, 1, dim)
        end
    end
    ℒ = Directional(f, v)
    return RadialBasisOperator(ℒ, data, eval_points, basis; k=k)
end

function RadialBasisOperator(
    ℒ::Directional,
    data::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    Na = length(adjl)
    Nd = length(data)
    weights = spzeros(eltype(D), Na, Nd)
    return RadialBasisOperator(ℒ, weights, data, data, adjl, basis)
end

function RadialBasisOperator(
    ℒ::Directional,
    data::AbstractVector{TD},
    eval_points::AbstractVector{TE},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {TD,TE,T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, eval_points, k)
    Na = length(adjl)
    Nd = length(data)
    weights = spzeros(eltype(TD), Na, Nd)
    return RadialBasisOperator(ℒ, weights, data, eval_points, adjl, basis)
end

function update_weights!(op::RadialBasisOperator{<:Directional})
    v = op.ℒ.v
    N = length(first(op.data))
    @assert length(v) == N || length(v) == size(op)[1] "wrong size for v"
    if length(v) == N
        op.weights .= mapreduce(+, enumerate(op.ℒ.ℒ)) do (i, ℒ)
            _build_weights(ℒ, op) * v[i]
        end
    else
        vv = ntuple(i -> getindex.(v, i), N)
        op.weights .= mapreduce(+, enumerate(op.ℒ.ℒ)) do (i, ℒ)
            Diagonal(vv[i]) * _build_weights(ℒ, op)
        end
    end
    validate_cache(op)
    return nothing
end

Base.size(op::RadialBasisOperator{<:Directional}) = size(op.weights)

# pretty printing
print_op(op::Directional) = "Directional Gradient (∇f⋅v)"
