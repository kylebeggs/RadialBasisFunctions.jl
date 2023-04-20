function ∂(m::MonomialBasis, order::T, dim::T) where {T<:Int}
    e = multiexponents(m.n + 1, m.deg)
    e = map(i -> getindex.(e, i), 1:(m.n))
    c = ones(T, length(e[dim]))
    ∂exponents!(e, c, order, dim)
    ids = map(j -> map(i -> findall(x -> x >= i, e[j]), 1:(m.deg)), 1:(m.n))
    return build_monomial_basis(m.n, m.deg, ids, c)
end

function ∇(m::MonomialBasis, dim::T) where {T<:Int}
    e = multiexponents(m.n + 1, m.deg)
    e = map(i -> getindex.(e, i), 1:(m.n))
    c = zeros(T, length(e[dim]))
    for i in eachindex(e[dim])
        if e[dim][i] > 0
            c[i] = e[dim][i]
            e[dim][i] -= 1
        end
    end
    ids = map(j -> map(i -> findall(x -> x >= i, e[j]), 1:(m.deg)), 1:(m.n))
    return build_monomial_basis(m.n, m.deg, ids, c)
end

∇(m::MonomialBasis) = ∇ℒ(x) = ForwardDiff.jacobian(m, x)

function ∇²(m::MonomialBasis)
    return function ∇²ℒ(x)
        return sum([ForwardDiff.jacobian(∂(m, dim), x)[:, dim] for dim in eachindex(x)])
    end
end

function build_monomial_basis(n::T, deg::T, ids, c::Vector{T}) where {T<:Int}
    function basis(x::AbstractVector{T}) where {T}
        b = ones(T, binomial(n + deg, n))
        # TODO flatten loop
        @views @inbounds for i in eachindex(ids), k in eachindex(ids[i])
            b[ids[i][k]] *= x[i]
        end
        return c .* b
    end
    return basis
end

function ∂exponents!(e, c, order, dim)
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
