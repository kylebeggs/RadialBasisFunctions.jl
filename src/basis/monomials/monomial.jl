"""
    struct MonomialBasis{T<:Int,B<:Function}

Multivariate Monomial basis.
n ∈ N: length of array, i.e., x ∈ Rⁿ
deg ∈ N: degree
"""
struct MonomialBasis{T<:Int,B<:Function}
    n::T
    deg::T
    basis::B
    function MonomialBasis(n::T, deg::T) where {T<:Int}
        basis = build_monomial_basis(n, deg)
        return new{T,typeof(basis)}(n, deg, basis)
    end
end

function (m::MonomialBasis)(x::AbstractVector{T}) where {T}
    b = ones(T, binomial(m.n + m.deg, m.n))
    m.basis(b, x)
    return b
end
(m::MonomialBasis)(b, x) = m.basis(b, x)

function Base.show(io::IO, pb::MonomialBasis)
    return print(io, "MonomialBasis, deg=$(pb.deg) in $(pb.n) dimensions")
end

function build_monomial_basis(n::T, deg::T) where {T<:Int}
    e = multiexponents(n + 1, deg)
    e = map(i -> getindex.(e, i), 1:n)
    ids = map(j -> map(i -> findall(x -> x >= i, e[j]), 1:deg), 1:n)
    return build_monomial_basis(ids)
end

function build_monomial_basis(ids::Vector{Vector{Vector{T}}}) where {T<:Int}
    function basis(b::AbstractVector{B}, x::AbstractVector) where {B}
        b .= one(B)
        # TODO flatten loop - why does it allocate here
        @views @inbounds for i in eachindex(ids), k in eachindex(ids[i])
            b[ids[i][k]] *= x[i]
        end
        return nothing
    end
    return basis
end
