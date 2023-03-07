struct Partial{S<:AbstractMatrix,V<:AbstractVector,T<:Int,B<:AbstractRadialBasis} <:
       ScalarValuedOperator
    ℒ::Function
    order::Int
    dim::Int
    weights::S
    data::Vector{V}
    adjl::Vector{Vector{T}}
    basis::B
    scales::Any
end

function Partial(
    ℒ::Function,
    order::Int,
    data::Vector,
    dim::Int,
    basis::B=PHS(3, 1);
    k::Int=autoselect_k(data, basis),
) where {B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    scales = find_scales(data, adjl)
    weights = build_weightmx(ℒ, data, adjl, basis, scales)
    return Partial(ℒ, order, dim, weights, data, adjl, basis, scales)
end

function Partial(
    data::Vector, dim::Int, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    ∂ℒ(x) = ∂(x, dim)
    return Partial(∂ℒ, 1, data, dim, basis; k=k)
end

function Partial²(
    data::Vector, dim::Int, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    ∂²ℒ(x) = ∂²(x, dim)
    return Partial(∂²ℒ, 2, data, dim, basis; k=k)
end

struct Laplacian{S<:AbstractMatrix,V,T<:Int,B<:AbstractRadialBasis} <: ScalarValuedOperator
    weights::S
    data::Vector{V}
    adjl::Vector{Vector{T}}
    basis::B
end

function Laplacian(
    data::Vector, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    n = length(data)
    weights = spzeros(n, n)
    for dim in eachindex(first(data))
        weights += build_weightmx(x -> ∂²(x, dim), data, adjl, basis)
    end
    return Laplacian(weights, data, adjl, basis)
end

# evaluate
(op::ScalarValuedOperator)(x::AbstractVecOrMat) = op.weights * x
(op::ScalarValuedOperator)(x::AbstractVector{<:AbstractVector}) = (op.weights,) .* x
function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::ScalarValuedOperator, x::AbstractVecOrMat
)
    return mul!(y, op.weights, x)
end
function LinearAlgebra.mul!(
    y::AbstractVecOrMat, op::ScalarValuedOperator, x::AbstractVecOrMat, α, β
)
    return mul!(y, op.weights, x, α, β)
end

# pretty printing
function Base.show(io::IO, op::Partial)
    println(io, "∂ⁿ/∂xᵢ (n = $(op.order), i = $(op.dim))")
    println(io, "  Data type: ", typeof(first(op.data)))
    println(io, "  Number of points: ", length(op.data))
    println(io, "  Dimensions: ", length(first(op.data)))
    println(io, "  Stencil size: ", length(first(op.adjl)))
    return println(io, "  Basis: ", op.basis, " with degree $(op.basis.deg) polynomial")
end

function Base.show(io::IO, op::Laplacian)
    println(io, "Laplacian, ∇²")
    println(io, "  Data type: ", typeof(first(op.data)))
    println(io, "  Number of points: ", length(op.data))
    println(io, "  Dimensions: ", length(first(op.data)))
    println(io, "  Stencil size: ", length(first(op.adjl)))
    return println(io, "  Basis: ", op.basis, " with degree $(op.basis.deg) polynomial")
end
