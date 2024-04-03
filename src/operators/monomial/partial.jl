# 1 dimensional, derivative w.r.t 1st (only) dimension
∂_get_ops(::Val{1}, ::Val{0}, ::Val{1}) = ((s, _) -> zero(s),)
function ∂_get_ops(::Val{1}, ::Val{N}, ::Val{1}) where {N}
    return (
        (s, _) -> zero(s),
        (s, _) -> s,
        ntuple(i -> ((s, x) -> s * (i + 1) / i * x[1]), N - 1)...,
    )
end

# 2 dimensional, derivative w.r.t 1st dimension
∂_get_ops(::Val{2}, ::Val{0}, ::Val{1}) = ((s, _) -> zero(s),)
function ∂_get_ops(::Val{2}, ::Val{1}, ::Val{1})
    return ((s, _) -> zero(s), (s, _) -> s, (s, _) -> zero(s))
end
function ∂_get_ops(::Val{2}, ::Val{2}, ::Val{1})
    return (
        (s, _) -> zero(s),
        (s, _) -> s,
        (s, x) -> s * 2 * x[1],
        (s, x) -> zero(s),
        (s, x) -> s * x[2] / (2 * x[1]),
        (s, x) -> zero(s),
    )
end

# 2 dimensional, derivative w.r.t 1st dimension
∂_get_ops(::Val{2}, ::Val{0}, ::Val{2}) = ((s, _) -> s,)
function ∂_get_ops(::Val{2}, ::Val{1}, ::Val{2})
    return ((s, _) -> zero(s), (s, _) -> zero(s), (s, _) -> s)
end
function ∂_get_ops(::Val{2}, ::Val{2}, ::Val{2})
    return (
        (s, _) -> zero(s),
        (s, _) -> zero(s),
        (s, _) -> zero(s),
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s / x[1] * 2 * x[2],
    )
end

## TODO: Implement the rest of the derivative operators
function ∂_get_ops(::Val{3}, ::Val{0})
    return ((s, _) -> s,)
end
function ∂_get_ops(::Val{3}, ::Val{1})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[3] / x[2],
    )
end
function ∂_get_ops(::Val{3}, ::Val{2})
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

∂²_get_ops(::Val{1}, ::Val{0}) = ((s, _) -> s,)
function ∂²_get_ops(::Val{1}, ::Val{N}) where {N}
    return ((s, _) -> s, ntuple(i -> ((s, x) -> s * x[1]), N)...)
end

∂²_get_ops(::Val{2}, ::Val{0}) = ((s, _) -> s,)
function ∂²_get_ops(::Val{2}, ::Val{1})
    return ((s, _) -> s, (s, x) -> s * x[1], (s, x) -> s * x[2] / x[1])
end
function ∂²_get_ops(::Val{2}, ::Val{2})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[1],
        (s, x) -> s / (x[1] * x[1]) * x[2],
        (s, x) -> s * x[1],
        (s, x) -> s / x[1] * x[2],
    )
end

function ∂²_get_ops(::Val{3}, ::Val{0})
    return ((s, _) -> s,)
end
function ∂²_get_ops(::Val{3}, ::Val{1})
    return (
        (s, _) -> s,
        (s, x) -> s * x[1],
        (s, x) -> s * x[2] / x[1],
        (s, x) -> s * x[3] / x[2],
    )
end
function ∂²_get_ops(::Val{3}, ::Val{2})
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
