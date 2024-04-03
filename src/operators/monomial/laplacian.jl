∇²_get_ops(::Val{1}, ::Val{0}) = ((s, _) -> s,)
function ∇²_get_ops(::Val{1}, ::Val{N}) where {N}
    return ((s, _) -> s, ntuple(i -> ((s, x) -> s * x[1]), N)...)
end

∇²_get_ops(::Val{2}, ::Val{0}) = ((s, _) -> s,)
function ∇²_get_ops(::Val{2}, ::Val{1})
    return ((s, _) -> s, (s, x) -> s * x[1], (s, x) -> s * x[2] / x[1])
end
function ∇²_get_ops(::Val{2}, ::Val{2})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[1],
        (s, x) -> s / (x[1] * x[1]) * x[2],
        (s, x) -> s * x[1],
        (s, x) -> s / x[1] * x[2],
    )
end

function ∇²_get_ops(::Val{3}, ::Val{0})
    return ((s, _) -> s,)
end
function ∇²_get_ops(::Val{3}, ::Val{1})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[3] / x[2],
    )
end
function ∇²_get_ops(::Val{3}, ::Val{2})
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
