struct Monomial{E,C}
    exponents::E
    coeffs::C
end

function ∂(mb::MonomialBasis, order::T, dim::T) where {T<:Int}
    #if mb.n <= 3 && mb.deg <= 2
    #    if order == 1
    #        return ∂(mb, dim)
    #    elseif order == 2
    #        return ∂²(mb, dim)
    #    else
    #        throw(
    #            ArgumentError(
    #                "order > 2 not currently supported with polynomial augmentation."
    #            ),
    #        )
    #    end
    #else
    me = ∂exponents(mb, order, dim)
    ids = monomial_recursive_list(mb, me)
    basis = build_monomial_basis(ids, me.coeffs)
    return basis
    #end
end

function ∂auto!(db, f, x, shadow)
    b = similar(db)
    Enzyme.autodiff_deferred(
        Forward, f, DuplicatedNoNeed(b, db), Duplicated(x, convert(typeof(x), shadow))
    )
    return nothing
end

function ∂(m::MonomialBasis, dim::Int)
    shadow = zeros(Int, m.n)
    shadow[dim] = 1
    function ∂ℒ!(db, x)
        ∂auto!(db, m, x, shadow)
        return nothing
    end
    return ∂ℒ!
end

function ∂²(m::MonomialBasis, dim::Int)
    shadow = zeros(Int, m.n)
    shadow[dim] = 1
    dm!(db, x) = ∂auto!(db, m, x, shadow)
    function ∂²ℒ!(db, x)
        ∂auto!(db, dm!, x, shadow)
        return nothing
    end
    return ∂²ℒ!
end

∇(m::MonomialBasis) = ∇ℒ(x) = ForwardDiff.jacobian(m, x)

function ∇²(m::MonomialBasis)
    return function ∇²ℒ(x)
        return sum([ForwardDiff.jacobian(∂(m, dim), x)[:, dim] for dim in eachindex(x)])
    end
end

function build_monomial_basis(ids::Vector{Vector{Vector{T}}}, c::Vector{T}) where {T<:Int}
    function basis(db::AbstractVector{B}, x::AbstractVector) where {B}
        db .= one(B)
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
