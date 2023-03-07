struct Gradient{S<:AbstractMatrix,V<:AbstractVector,T<:Int,B<:AbstractRadialBasis} <:
       VectorValuedOperator
    weights::Vector{S}
    data::Vector{V}
    adjl::Vector{Vector{T}}
    basis::B
end

function Gradient(
    data::Vector, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    dims = 1:length(first(data))
    weights = [build_weightmx(x -> ∂(x, dim), data, adjl, basis) for dim in dims]
    return Gradient(weights, data, adjl, basis)
end

# TODO - finish hessian
function hessian(
    data::Vector, dim::Int, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    error("not implemented yet")
    return GenericRadialBasisOperator(x -> ∂²(x, dim), data, basis; k=k)
end

# evaluate
# TODO thread this
function (op::VectorValuedOperator)(x::AbstractVector{<:Real})
    return [op.weights[i] * x for i in eachindex(op.weights)]
end

(op::VectorValuedOperator)(x::AbstractVector{<:AbstractVector}) = op.weights .* x

function LinearAlgebra.:⋅(op::VectorValuedOperator, x::AbstractVector)
    return sum(op(x))
end

# TODO thread this
function LinearAlgebra.mul!(
    y::AbstractVector{<:Real}, op::VectorValuedOperator, x::AbstractVector{<:Real}
)
    for i in eachindex(op.weights)
        mul!(y[i], op.weights[i], x)
    end
end

# pretty printing
function Base.show(io::IO, op::Gradient)
    println(io, "Gradient")
    println(io, "  Data type: ", typeof(first(op.data)))
    println(io, "  Number of points: ", length(op.data))
    println(io, "  Dimensions: ", length(first(op.data)))
    println(io, "  Stencil size: ", length(first(op.adjl)))
    return println(io, "  Basis: ", op.basis, " with degree $(op.basis.deg) polynomial")
end
