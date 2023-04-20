abstract type VectorValuedOperator <: AbstractRadialBasisOperator end

struct Gradient{L<:Function} <: VectorValuedOperator
    ℒ::L
end

function gradient(
    data, basis::B=PHS(3; poly_deg=2); k::T=autoselect_k(data, basis)
) where {T<:Int,B<:AbstractRadialBasis}
    adjl = find_neighbors(data, k)
    weights = ntuple(
        dim -> _build_weightmx(x -> ∂(x, Val{1}(), dim), data, adjl, basis),
        length(first(data)),
    )
    return RadialBasisOperator(ℒ, weights, data, adjl, basis)
end

(op::VectorValuedOperator)(x) = op.ℒ(x)

# evaluate
function (op::RadialBasisOperator{VectorValuedOperator})(x::AbstractVector)
    return map(w -> w * x, op.weights)
end
function LinearAlgebra.:⋅(op::RadialBasisOperator{VectorValuedOperator}, x::AbstractVector)
    return sum(op(x))
end

# TODO
function LinearAlgebra.mul!(
    y::AbstractVector{<:Real}, op::VectorValuedOperator, x::AbstractVector
)
    for i in eachindex(op.weights)
        mul!(y[i], op.weights[i], x)
    end
end

# pretty printing
function Base.show(io::IO, op::Gradient)
    println(io, "Gradient")
    println(io, "└─Data type: ", typeof(first(op.data)))
    println(io, "└─Number of points: ", length(op.data))
    println(io, "└─Dimensions: ", length(first(op.data)))
    println(io, "└─Stencil size: ", length(first(op.adjl)))
    print(io, "└─RBF: ", op.basis)
    return nothing
end
