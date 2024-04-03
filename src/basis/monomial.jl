"""
    struct MonomialBasis{Dim,Deg}

`Dim` dimensional monomial basis of order `Deg`.
"""
struct MonomialBasis{Dim,Deg,Ops}
    ops::Ops
end
MonomialBasis{Dim,Deg}(ops::Ops) where {Dim,Deg,Ops} = MonomialBasis{Dim,Deg,Ops}(ops)
function MonomialBasis(dim::T, deg::T) where {T<:Int}
    ops = _get_ops(Val(dim), Val(deg))
    return MonomialBasis{dim,deg,typeof(ops)}(ops)
end
(m::MonomialBasis)(x) = MonomialBasisIterator(m, x)

struct MonomialBasisIterator{T,Dim,Deg,X}
    basis::MonomialBasis{Dim,Deg}
    x::X
    function MonomialBasisIterator(mb::MonomialBasis{Dim,Deg}, x::X) where {Dim,Deg,X}
        !_check_input(x) && throw(ArgumentError("Input must be a tuple of the same type"))
        return new{eltype(x),Dim,Deg,X}(mb, x)
    end
end

_check_input(::AbstractVector) = true
_check_input(x::Tuple) = all(a -> typeof(x[1]) == typeof(a), x)

mutable struct MonomialState{TV,TI}
    val::TV
    index::TI
end

function Base.iterate(m::MonomialBasisIterator{T}, state=MonomialState(one(T), 1)) where {T}
    if state.index > length(m.basis.ops)
        return nothing
    else
        val = m.basis.ops[state.index](state.val, m.x)
        state.val = iszero(val) ? state.val : val
        state.index += 1
        return (val, state)
    end
end
Base.length(m::MonomialBasisIterator) = length(m.basis.ops)
Base.eltype(::Type{<:MonomialBasisIterator{T}}) where {T} = T

function Base.show(io::IO, ::MonomialBasis{Dim,Deg}) where {Dim,Deg}
    return print(io, "MonomialBasis, deg=$(Deg) in $(Dim) dimensions")
end

_get_ops(::Val{1}, ::Val{0}) = ((s, _) -> s,)
function _get_ops(::Val{1}, ::Val{N}) where {N}
    return ((s, _) -> s, ntuple(i -> ((s, x) -> s * x[1]), N)...)
end

_get_ops(::Val{2}, ::Val{0}) = ((s, _) -> s,)
_get_ops(::Val{2}, ::Val{1}) = ((s, _) -> s, (s, x) -> s * x[1], (s, x) -> s * x[2] / x[1])
function _get_ops(::Val{2}, ::Val{2})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[1],
        (s, x) -> s / (x[1] * x[1]) * x[2],
        (s, x) -> s * x[1],
        (s, x) -> s / x[1] * x[2],
    )
end

function _get_ops(::Val{3}, ::Val{0})
    return ((s, _) -> s,)
end
function _get_ops(::Val{3}, ::Val{1})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[3] / x[2],
    )
end
function _get_ops(::Val{3}, ::Val{2})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[3] / x[2],
        (s, x) -> s * x[1] * x[1] / x[3],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[1] * x[3] / (x[2] * x[2]),
        (s, x) -> s * x[3] / x[1],
        (s, x) -> s * x[2] / x[3],
    )
end
