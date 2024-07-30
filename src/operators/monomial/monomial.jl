include("partial.jl")
∂(mb::MonomialBasis, differentiation_dim::Int) = ∂(mb, Val(differentiation_dim))
∂²(mb::MonomialBasis, differentiation_dim::Int) = ∂²(mb, Val(differentiation_dim))

#function ∇(m::MonomialBasis{Dim}) where {Dim}
#    basis(x) = ntuple(dim -> ∂(m, dim)(x), Dim)
#    return ℒMonomialBasis(m,basis)
#end

#function ∇²(m::MonomialBasis{Dim}) where {Dim}
#    basis(x) = sum(dim -> ∂²(m, dim)(x), 1:Dim)
#    return ℒMonomialBasis(m,basis)
#end

#########################################################################

∂(basis::MonomialBasis, ::Val{N}) where {N} = _∂(basis, 1, Val(N))
∂²(basis::MonomialBasis, ::Val{N}) where {N} = _∂(basis, 2, Val(N))

function _∂(m::MonomialBasis{Dim,Deg}, order::Int, ::Val{N}) where {Dim,Deg,N}
    me = ∂exponents(m, order, N)
    ids = monomial_recursive_list(m, me)
    basis! = build_monomial_basis(ids, me.coeffs)
    return ℒMonomialBasis(Dim, Deg, basis!)
end

function ∇²(m::MonomialBasis{Dim,Deg}) where {Dim,Deg}
    ∂² = ntuple(dim -> ∂(m, 2, dim), Dim)
    function basis!(b, x)
        cache = ones(size(b))
        b .= 0
        for ∂²! in ∂²
            # use mapreduce here instead?
            ∂²!(cache, x)
            b .+= cache
        end
        return nothing
    end
    return ℒMonomialBasis(Dim, Deg, basis!)
end

struct Monomial{E,C}
    exponents::E
    coeffs::C
end

function build_monomial_basis(ids::Vector{Vector{Vector{T}}}, c::Vector{T}) where {T<:Int}
    function basis!(db::AbstractVector{B}, x::AbstractVector) where {B}
        db .= 1
        # TODO flatten loop - why does it allocate here
        @views @inbounds for i in eachindex(ids), j in eachindex(ids[i])
            db[ids[i][j]] *= x[i]
        end
        db .*= c
        return nothing
    end
    return basis!
end

function ∂exponents(::MonomialBasis{Dim,Deg}, order::T, dim::T) where {Dim,Deg,T<:Int}
    ex = collect(Vector{Int}, multiexponents(Dim + 1, Deg))
    N = binomial(Dim + Deg, Deg)
    e = [zeros(Int, N) for _ in 1:Dim]
    for i in 1:(Dim), j in 1:N
        e[i][j] = ex[j][i]
    end
    c = ones(T, length(e[dim]))
    ∂exponents!(e, c, order, dim)
    return Monomial(e, c)
end

function ∂exponents!(e, c, order::Int, dim::Int)
    order == 0 && return nothing
    for i in eachindex(e[dim])
        e[dim][i] == 0 && (c[i] = 0)
        if e[dim][i] > 0
            c[i] *= e[dim][i]
            e[dim][i] -= 1
        end
    end
    return ∂exponents!(e, c, order - 1, dim)
end

function monomial_recursive_list(::MonomialBasis{Dim,Deg}, me::Monomial) where {Dim,Deg}
    return Vector{Vector{Int}}[
        Vector{Int}[findall(x -> x >= i, me.exponents[j]) for i in 1:(Deg)] for j in 1:(Dim)
    ]
end

function monomial_recursive_list(m::MonomialBasis, me::Vector{<:Monomial})
    return [monomial_recursive_list(m, me[i]) for i in eachindex(me)]
end
