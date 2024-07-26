"""
    struct MonomialBasis{Dim,Deg} <: AbstractBasis

`Dim` dimensional monomial basis of order `Deg`.
"""
struct MonomialBasis{Dim,Deg} <: AbstractBasis
    function MonomialBasis(dim::T, deg::T) where {T<:Int}
        if deg < 0
            throw(ArgumentError("Monomial basis must have non-negative degree"))
        end
        return new{dim,deg}()
    end
end

(m::MonomialBasis{1,0})(x) = one(_get_underlying_type(x))
(m::MonomialBasis{1,0})(x::AbstractVector) = SVector{1}(one(_get_underlying_type(x)))
function (m::MonomialBasis{1,N})(x) where {N}
    return SVector{N + 1}(one(x), ntuple(n -> x^n, N)...)
end
function (m::MonomialBasis{1,N})(x::AbstractVector) where {N}
    return SVector{N + 1}(one(eltype(x)), ntuple(n -> x[1]^n, N)...)
end

(m::MonomialBasis{2,0})(x) = one(eltype(x))
(m::MonomialBasis{2,1})(x) = SVector{3}(one(eltype(x)), x[1], x[2])
function (m::MonomialBasis{2,2})(x)
    return SVector{6}(one(eltype(x)), x[1], x[2], x[1] * x[2], x[1] * x[1], x[2] * x[2])
end

(m::MonomialBasis{3,0})(x) = one(eltype(x))
(m::MonomialBasis{3,1})(x) = SVector{4}(one(eltype(x)), x[1], x[2], x[3])
function (m::MonomialBasis{3,2})(x)
    return SVector{10}(
        one(eltype(x)),
        x[1],
        x[2],
        x[3],
        x[1] * x[2],
        x[1] * x[3],
        x[2] * x[3],
        x[1] * x[1] * x[1],
        x[2] * x[2] * x[2],
        x[3] * x[3] * x[3],
    )
end

function Base.show(io::IO, ::MonomialBasis{Dim,Deg}) where {Dim,Deg}
    return print(io, "MonomialBasis of degree $(Deg) in $(Dim) dimensions")
end
