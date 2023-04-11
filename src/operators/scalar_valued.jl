abstract type ScalarValuedOperator <: AbstractRadialBasisOperator end

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
    f = let o = Val{order}, dim = dim
        x -> ∂(x, o, dim)
    end
    ℒ = Partial(f, order, dim)
    adjl = find_neighbors(data, k)
    N = length(data)
    weights = spzeros(N, N)
    return RadialBasisOperator(ℒ, weights, data, adjl, basis)
end

"""
    Laplacian <: ScalarValuedOperator

Builds an operator for the sum of the second derivatives w.r.t. each independent variable.
"""
struct Laplacian{L<:Function} <: ScalarValuedOperator
    ℒ::L
end

# convienience constructors
function laplacian(
    data::AbstractVector{D}, basis::B=PHS(3; poly_deg=2); k::T=autoselect_k(data, basis)
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    ℒ = Laplacian(∇²)
    adjl = find_neighbors(data, k)
    N = length(data)
    weights = spzeros(N, N)
    return RadialBasisOperator(ℒ, weights, data, adjl, basis)
end

(op::ScalarValuedOperator)(x) = op.ℒ(x)

# pretty printing
print_op(op::Partial) = "∂ⁿ/∂xᵢ (n = $(op.order), i = $(op.dim))"
print_op(op::Laplacian) = "Laplacian (∇² or Δ)"
