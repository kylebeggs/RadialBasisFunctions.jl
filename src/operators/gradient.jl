"""
    Gradient <: VectorValuedOperator

Builds an operator for the gradient of a function.
"""
struct Gradient{L<:NTuple} <: VectorValuedOperator
    ℒ::L
end

# convienience constructors
"""
    function gradient(data, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the gradient, `Gradient`.
"""
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

"""
    function gradient(data, eval_points, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the gradient, `Gradient`. The resulting operator will only evaluate at `eval_points`.
"""
function gradient(
    data::AbstractVector{D},
    eval_points::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,B<:AbstractRadialBasis,T<:Int}
    f = ntuple(length(first(data))) do dim
        return let dim = dim
            x -> ∂(x, 1, dim)
        end
    end
    ℒ = Gradient(f)
    return RadialBasisOperator(ℒ, data, eval_points, basis; k=k)
end

function RadialBasisOperator(
    ℒ::Gradient,
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

function RadialBasisOperator(
    ℒ::Gradient,
    data::AbstractVector{D},
    eval_points::AbstractVector{D},
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    TD = eltype(D)
    adjl = find_neighbors(data, eval_points, k)
    Na = length(adjl)
    Nd = length(data)
    weights = ntuple(_ -> _allocate_weights(TD, Na, Nd, k), length(ℒ.ℒ))
    return RadialBasisOperator(ℒ, weights, data, eval_points, adjl, basis)
end

Base.size(op::RadialBasisOperator{<:Gradient}) = size(first(op.weights))

# pretty printing
print_op(op::Gradient) = "Gradient (∇f)"
