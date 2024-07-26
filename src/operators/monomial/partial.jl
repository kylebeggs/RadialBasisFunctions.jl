function ∂(::MonomialBasis{Dim,0}, ::Val{N}) where {Dim,N}
    basis(x) = zero(_get_underlying_type(x))
    return ℒMonomial(basis)
end

_∂(x, n) = n * x^(n - 1)
_∂(x::Union{AbstractArray,Tuple}, n) = n * x[1]^(n - 1)
function ∂(::MonomialBasis{1,Deg}, ::Val{N}) where {Deg,N}
    function basis(x)
        x_type = _get_underlying_type(x)
        return SVector{Deg + 1}(zero(x_type), ntuple(Base.Fix1(_∂, x), Deg)...)
    end
    return ℒMonomial(basis)
end
function ∂(::MonomialBasis{1,0}, ::Val{N}) where {N}
    basis(x) = zero(_get_underlying_type(x))
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{2,1}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{3}(zero(T), one(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{2,1}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{3}(zero(T), zero(T), one(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{2,2}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{6}(zero(T), one(T), zero(T), x[2], 2 * x[1], zero(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{2,2}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{6}(zero(T), zero(T), one(T), x[1], zero(T), 2 * x[2])
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,1}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), one(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,1}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), zero(T), one(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,1}, ::Val{3})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), zero(T), zero(T), one(T))
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,2}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            one(T),
            zero(T),
            zero(T),
            x[2],
            x[3],
            zero(T),
            3 * x[1] * x[1],
            zero(T),
            zero(T),
        )
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,2}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            zero(T),
            one(T),
            zero(T),
            x[1],
            zero(T),
            x[3],
            zero(T),
            3 * x[2] * x[2],
            zero(T),
        )
    end
    return ℒMonomial(basis)
end

function ∂(::MonomialBasis{3,2}, ::Val{3})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            zero(T),
            zero(T),
            one(T),
            zero(T),
            x[1],
            x[2],
            zero(T),
            zero(T),
            3 * x[3] * x[3],
        )
    end
    return ℒMonomial(basis)
end

## ∂²

∂²(m::MonomialBasis{Dim,0}, n::Val{N}) where {Dim,N} = ∂(m, n)
∂²(m::MonomialBasis{1,0}, n::Val{N}) where {N} = ∂(m, n)

_∂²(x, n) = (n - 1) * n * x^(n - 2)
_∂²(x::Union{AbstractArray,Tuple}, n) = (n - 1) * n * x[1]^(n - 2)
function ∂²(::MonomialBasis{1,Deg}, ::Val{N}) where {Deg,N}
    function basis(x)
        x_type = _get_underlying_type(x)
        return SVector{Deg + 1}(zero(x_type), ntuple(Base.Fix1(_∂², x), Deg)...)
    end
    return ℒMonomial(basis)
end
function ∂²(::MonomialBasis{1,1}, ::Val{N}) where {N}
    function basis(x)
        x_type = _get_underlying_type(x)
        return SVector{2}(zero(x_type), zero(x_type))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{2,1}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{3}(zero(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{2,1}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{3}(zero(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{2,2}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{6}(zero(T), zero(T), zero(T), zero(T), 2, zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{2,2}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{6}(zero(T), zero(T), zero(T), zero(T), zero(T), 2)
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,1}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), zero(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,1}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), zero(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,1}, ::Val{3})
    function basis(x)
        T = eltype(x)
        return SVector{4}(zero(T), zero(T), zero(T), zero(T))
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,2}, ::Val{1})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            6 * x[1],
            zero(T),
            zero(T),
        )
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,2}, ::Val{2})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            6 * x[2],
            zero(T),
        )
    end
    return ℒMonomial(basis)
end

function ∂²(::MonomialBasis{3,2}, ::Val{3})
    function basis(x)
        T = eltype(x)
        return SVector{10}(
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            zero(T),
            6 * x[3],
        )
    end
    return ℒMonomial(basis)
end
