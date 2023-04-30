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
        #if n <= 3 && deg <= 2
        #    basis = build_monomial_basis(Val{n}(), Val{deg}())
        #else
        #    basis = build_monomial_basis(n, deg)
        #end

        # TODO fix Enzyme hanging
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

function build_monomial_basis(::Any, ::Val{0})
    function basis(b::AbstractVector, ::AbstractVector{T}) where {T}
        b[begin] = one(T)
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{1}, ::Val{1})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{1}, ::Val{2})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        b[begin + 2] = x[1]^2
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{2}, ::Val{1})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        b[begin + 2] = x[2]
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{2}, ::Val{2})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        b[begin + 2] = x[2]
        b[begin + 3] = x[1]^2
        b[begin + 4] = x[2]^2
        b[begin + 5] = x[1] * x[2]
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{3}, ::Val{1})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        b[begin + 2] = x[2]
        b[begin + 3] = x[3]
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{3}, ::Val{2})
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        b[begin] = one(T)
        b[begin + 1] = x[1]
        b[begin + 2] = x[2]
        b[begin + 3] = x[3]
        b[begin + 4] = x[1]^2
        b[begin + 5] = x[2]^2
        b[begin + 6] = x[3]^2
        b[begin + 7] = x[1] * x[2]
        b[begin + 8] = x[1] * x[3]
        b[begin + 9] = x[2] * x[3]
        return nothing
    end
    return basis
end

function build_monomial_basis(::Val{2}, deg::T) where {T<:Int}
    # pascals triangle
    function basis(b::AbstractVector, x::AbstractVector{T}) where {T}
        offset = 0
        for line in 1:(d + 1)
            for k in 1:(line - 1)
                for i in 1:k
                    b[i + offset] *= x[1]
                end
                for i in (line - k + 1):line
                    b[i + offset] *= x[2]
                end
            end
            offset += line
        end
        return nothing
    end
    return basis
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
