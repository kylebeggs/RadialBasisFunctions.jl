"""
    Directional <: VectorValuedOperator

Operator for the directional derivative, or the inner product of the gradient and a direction vector.
"""
struct Directional{L<:NTuple,T} <: VectorValuedOperator
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

function RadialBasisOperator(
    ℒ::Directional,
    data::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    TD = eltype(D)
    adjl = find_neighbors(data, k)
    N = length(adjl)
    weights = ntuple(_ -> _allocate_weights(TD, N, N, k), length(ℒ.ℒ))
    return RadialBasisOperator(ℒ, weights, data, data, adjl, basis)
end

function _update_weights!(
    op::RadialBasisOperator{<:Directional}, weights::NTuple{N,AbstractMatrix}
) where {N}
    v = op.ℒ.v
    @assert length(v) == N || length(v) == size(op)[2] "wrong size for v"
    if length(v) == N
        for (i, ℒ) in enumerate(op.ℒ.ℒ)
            weights[i] .= _build_weights(ℒ, op) * v[i]
        end
    else
        vv = ntuple(i -> getindex.(v, i), N)
        for (i, ℒ) in enumerate(op.ℒ.ℒ)
            weights[i] .= Diagonal(vv[i]) * _build_weights(ℒ, op)
        end
    end
    validate_cache(op)
    return nothing
end

Base.size(op::RadialBasisOperator{<:Directional}) = size(first(op.weights))

# pretty printing
print_op(op::Directional) = "Directional Gradient (∇f⋅v)"
