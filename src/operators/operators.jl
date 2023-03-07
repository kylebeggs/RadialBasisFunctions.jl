"""
    struct GenericRadialBasisOperator{S <: AbstractMatrix, V, T <: Int}

Operator of data using a radial basis with potential polynomial augmentation.
"""
struct GenericRadialBasisOperator{
    F,S<:AbstractMatrix,V<:AbstractVector,T<:Int,B<:AbstractRadialBasis
} <: AbstractRadialBasisOperator
    ℒ::Function
    weights::S
    data::Vector{V}
    adjl::Vector{Vector{T}}
    basis::B
end

function GenericRadialBasisOperator(
    ℒ::Function, data::Vector, basis::B=PHS(3, 2); k::Int=autoselect_k(data, basis)
) where {B<:AbstractRadialBasis}
    # TODO
    #check_num_params(ℒ, basis)
    adjl = find_neighbors(data, k)
    weights = build_weightmx(ℒ, data, adjl, basis)
    return GenericRadialBasisOperator(ℒ, weights, data, adjl, basis)
end
Base.getindex(op::GenericRadialBasisOperator, i) = nonzeros(op.weights[i, :])

# include built-in operators
include("scalar_valued.jl")
include("vector_valued.jl")

# pretty printing
function Base.show(io::IO, op::GenericRadialBasisOperator)
    println(io, "GenericRadialBasisOperator")
    println(io, "  Data type: ", typeof(first(op.data)))
    println(io, "  Number of points: ", length(op.data))
    println(io, "  Dimensions: ", length(first(op.data)))
    println(io, "  Stencil size: ", length(first(op.adjl)))
    return println(
        io, "  Basis: ", typeof(op.basis), " with degree $(op.basis.deg) polynomial"
    )
end
