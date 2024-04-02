"""
    struct MonomialBasis{Dim,Deg}

Multivariate Monomial basis.
n ∈ N: length of array, i.e., x ∈ Rⁿ
deg ∈ N: degree
"""
struct MonomialBasis{Dim,Deg}
    MonomialBasis(dim::T, deg::T) where {T<:Int} = new{dim,deg}()
end
(m::MonomialBasis{Dim,Deg})(x) where {Dim,Deg} = MonomialBasisIterator(x, Deg)

struct MonomialBasisIterator{T,D,X,O}
    x::X
    ops::O
    function MonomialBasisIterator{D}(x::X, ops::O) where {D,X,O}
        !_check_input(x) && throw(ArgumentError("Input must be a tuple of the same type"))
        return new{eltype(x),D,X,O}(x, ops)
    end
end

function MonomialBasisIterator(x, deg::T) where {T<:Int}
    N = length(x)
    ops = _get_ops(Val(N), Val(deg))
    return MonomialBasisIterator{deg}(x, ops)
end

_check_input(::AbstractVector) = true
_check_input(x::Tuple) = all(a -> typeof(x[1]) == typeof(a), x)

function _get_ops(::Val{2}, ::Val{0})
    return ((s, _) -> s,)
end

function _get_ops(::Val{2}, ::Val{1})
    return ((s, _) -> s, (s, x) -> s * x[1], (s, x) -> s * x[2] / x[1])
end

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

mutable struct MonomialState{TV,TI}
    val::TV
    index::TI
end

function Base.iterate(m::MonomialBasisIterator{T}, state=MonomialState(one(T), 1)) where {T}
    if state.index > length(m.ops)
        return nothing
    else
        state.val = m.ops[state.index](state.val, m.x)
        state.index += 1
        return (state.val, state)
    end
end
Base.length(m::MonomialBasisIterator) = length(m.ops)
Base.eltype(::Type{<:MonomialBasisIterator{T}}) where {T} = T

# pretty printing
function Base.show(io::IO, pb::MonomialBasis)
    return print(io, "MonomialBasis, deg=$(pb.deg) in $(pb.n) dimensions")
end
