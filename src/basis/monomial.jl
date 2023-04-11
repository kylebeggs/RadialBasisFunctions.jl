"""
    struct MonomialBasis{N<:Int,D<:Int,B}

Multivariate Monomial basis.
n ∈ N: length of array, i.e., x ∈ Rⁿ
deg ∈ N: degree
"""
struct MonomialBasis{T<:Int,B<:Function}
    n::T
    deg::T
    basis::B
    function MonomialBasis(n::T, deg::T) where {T<:Int}
        if n <= 3 && deg <= 2
            basis = m_basis(Val{n}, Val{deg})
        else
            basis = m_basis(n, deg)
        end
        return new{T,typeof(basis)}(n, deg, basis)
    end
end

(m::MonomialBasis)(x) = m.basis(x)

∂(m::MonomialBasis, dim::Int) = ∂ℒ(x) = ForwardDiff.jacobian(m, x)[:, dim]
∇(m::MonomialBasis) = ∇ℒ(x) = ForwardDiff.jacobian(m, x)
∂²(m::MonomialBasis, dim::Int) = ∂²ℒ(x) = ForwardDiff.jacobian(∂(m, dim), x)[:, dim]
function ∇²(m::MonomialBasis)
    return function ∇²ℒ(x)
        return sum([ForwardDiff.jacobian(∂(m, dim), x)[:, dim] for dim in eachindex(x)])
    end
end

function Base.show(io::IO, pb::MonomialBasis)
    return print(io, "MonomialBasis, deg=$(pb.deg) in $(pb.n) dimensions")
end

function m_basis(n::T, deg::T) where {T<:Int}
    e = multiexponents(n + 1, deg)
    e = map(i -> getindex.(e, i), 1:n)
    idxs = map(j -> map(i -> findall(x -> x >= i, e[j]) .+ 1, 1:deg), 1:n)
    function basis(x::T) where {T}
        b = ones(eltype(T), binomial(n + deg, n))
        # TODO flatten loop
        @inbounds for i in eachindex(idxs), k in eachindex(idxs[i])
            b[idxs[i][k]] *= x[i]
        end
        return b
    end
    return basis
end

function m_basis(::Any, ::Type{Val{0}})
    return basis(::AbstractVector{T}) where {T} = one(T)
end

function m_basis(::Type{Val{1}}, ::Type{Val{1}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{2,T}(one(T), x[1])
    end
    return basis
end

function m_basis(::Type{Val{1}}, ::Type{Val{2}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{3,T}(one(T), x[1], x[1]^2)
    end
    return basis
end

function m_basis(::Type{Val{2}}, ::Type{Val{1}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{3,T}(one(T), x[1], x[2])
    end
    return basis
end

function m_basis(::Type{Val{2}}, ::Type{Val{2}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{6,T}(one(T), x[1], x[2], x[1]^2, x[2]^2, x[1] * x[2])
    end
    return basis
end

function m_basis(::Type{Val{3}}, ::Type{Val{1}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{4,T}(one(T), x[1], x[2], x[3])
    end
    return basis
end

function m_basis(::Type{Val{3}}, ::Type{Val{2}})
    function basis(x::AbstractVector{T}) where {T}
        return SVector{10,T}(
            one(T),
            x[1],
            x[2],
            x[3],
            x[1]^2,
            x[2]^2,
            x[3]^2,
            x[1] * x[2],
            x[2] * x[3],
            x[1] * x[3],
        )
    end
    return basis
end

function pascals_triangle(x::T, d::N) where {T,N<:Int}
    n = length(x)
    c = ones(eltype(T), binomial(n + d, n))
    offset = 0
    for line in 1:(d + 1)
        for k in 1:(line - 1)
            for i in 1:k
                c[i + offset] *= x[1]
            end
            for i in (line - k + 1):line
                c[i + offset] *= x[2]
            end
        end
        offset += line
    end
    return c
end

function pascals_tetrahedral(x::T, d::N) where {T,N<:Int}
    n = length(x)
    B = binomial(n + d, n)
    c = ones(eltype(T), B)
    for i in 1:(n - 2)
        for j in (n - 1):-1:0
            r!(c, x[i], i, i + j, d + 1, round(Int, B / 2))
        end
    end
    #for dd in 1:d
    #    build_level!(c, x, n, dd, binomial(n + dd - 1, n))
    #end
    return c
end

function build_level!(c, x, n, d, offset)
    for i in 1:(n - 2)
        for j in (n - 1):-1:0
            r!(c, x[i], i, i + j, d + 1, offset)
        end
    end
    return nothing
end

function r!(c, x, start, stop, k, offset)
    println("$start - $stop")
    for i in start:stop
        c[i + offset] *= x
    end
    if start != stop
        start += k
        stop += k - 1
        r!(c, x, start, stop, k, offset)
    else
        return nothing
    end
end
