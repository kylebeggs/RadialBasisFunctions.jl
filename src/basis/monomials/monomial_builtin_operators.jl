struct Monomial{E,C}
    exponents::E
    coeffs::C
end

function ∂(mb::MonomialBasis, order::T, dim::T) where {T<:Int}
    me = ∂exponents(mb, order, dim)
    ids = monomial_recursive_list(mb, me)
    basis = build_monomial_basis(ids, me.coeffs)
    return basis
end

∇(m::MonomialBasis) = ∇ℒ(x) = ForwardDiff.jacobian(m, x)

function ∇²(m::MonomialBasis)
    ∂² = ntuple(dim -> ∂(m, 2, dim), m.n)
    function ∇²ℒ(b, x)
        cache = ones(size(b))
        b .= 0
        for ∂²! in ∂²
            # use mapreduce here instead?
            ∂²!(cache, x)
            b .+= cache
        end
        return nothing
    end
    return ∇²ℒ
end

function build_monomial_basis(ids::Vector{Vector{Vector{T}}}, c::Vector{T}) where {T<:Int}
    function basis(db::AbstractVector{B}, x::AbstractVector) where {B}
        db .= 1
        # TODO flatten loop - why does it allocate here
        @views @inbounds for i in eachindex(ids), j in eachindex(ids[i])
            db[ids[i][j]] *= x[i]
        end
        db .*= c
        return nothing
    end
    return basis
end

function ∂exponents(mb::MonomialBasis, order::T, dim::T) where {T<:Int}
    ex = collect(Vector{Int}, multiexponents(mb.n + 1, mb.deg))
    N = binomial(mb.n + mb.deg, mb.deg)
    e = [zeros(Int, N) for _ in 1:(mb.n)]
    for i in 1:(mb.n), j in 1:N
        e[i][j] = ex[j][i]
    end
    c = ones(T, length(e[dim]))
    ∂exponents!(e, c, order, dim)
    return Monomial(e, c)
end

function ∂exponents!(e, c, order::T, dim::T) where {T<:Int}
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

function monomial_recursive_list(mb::MonomialBasis, me::Monomial)
    return Vector{Vector{Int}}[
        Vector{Int}[findall(x -> x >= i, me.exponents[j]) for i in 1:(mb.deg)] for
        j in 1:(mb.n)
    ]
end

function monomial_recursive_list(mb::MonomialBasis, me::Vector{<:Monomial})
    return [monomial_recursive_list(mb, me[i]) for i in eachindex(me)]
end
