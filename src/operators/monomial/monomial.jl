include("partial.jl")
∂(mb::MonomialBasis, differentiation_dim::Int) = ∂(mb, Val(differentiation_dim))
∂²(mb::MonomialBasis, differentiation_dim::Int) = ∂²(mb, Val(differentiation_dim))

function ∇(m::MonomialBasis{Dim}) where {Dim}
    basis(x) = ntuple(dim -> ∂(m, dim)(x), Dim)
    return ℒMonomial(basis)
end

function ∇²(m::MonomialBasis{Dim}) where {Dim}
    basis(x) = sum(dim -> ∂²(m, dim)(x), 1:Dim)
    return ℒMonomial(basis)
end
