########################################################################################
# Polyharmonic Spline

"""
    PHS

Polyharmonic spline , ϕ(r) = rⁿ. Allowable values of n are 1, 3, 5, or 7.
"""
struct PHS{N<:Int,D<:Int} <: AbstractRadialBasis
    n::N
    deg::D
    function PHS(n::N=3; deg::D=2) where {N<:Int,D<:Int}
        iseven(n) && throw(ArgumentError("n should be odd (n=$n)."))
        deg < -1 && throw(
            ArgumentError("Augmented polynomial degree must be at least 0 (constant).")
        )
        return new{N,D}(n, deg)
    end
end

# Dispatch on PHS order
(phs::PHS)(x, xᵢ) = ϕphs(Val{phs.n}, x, xᵢ)
∂(phs::PHS, dim::Int) = ∂phs(Val{phs.n}, dim)
∇(phs::PHS) = ∇phs(Val{phs.n})
∂²(phs::PHS, dim::Int) = ∂²phs(Val{phs.n}, dim)

# Linear polyharmonic spline, ϕ(r) = r.
ϕphs(::Type{Val{1}}, x, xᵢ) = euclidean(x, xᵢ)
∂phs(::Type{Val{1}}, dim::Int) = ∂ℒ(x, xᵢ) = (x[dim] - xᵢ[dim]) / euclidean(x, xᵢ)
∇phs(::Type{Val{1}}) = ∇ℒ(x, xᵢ) = (x .- xᵢ) / euclidean(x, xᵢ)
function ∂²phs(::Type{Val{1}}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return (-(x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ)) /
               (sqeuclidean(x, xᵢ)^1.5 + 1e-8)
    end
end

# Cubic polyharmonic spline, ϕ(r) = r³.
ϕphs(::Type{Val{3}}, x, xᵢ) = euclidean(x, xᵢ)^3
∂phs(::Type{Val{3}}, dim::Int) = ∂ℒ(x, xᵢ) = 3 * (x[dim] - xᵢ[dim]) * euclidean(x, xᵢ)
∇phs(::Type{Val{3}}) = ∇ℒ(x, xᵢ) = 3 * (x .- xᵢ) * euclidean(x, xᵢ)
function ∂²phs(::Type{Val{3}}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return 3 * (sqeuclidean(x, xᵢ) + (x[dim] - xᵢ[dim])^2) / (euclidean(x, xᵢ) + 1e-8)
    end
end

# Quintic polyharmonic spline, ϕ(r) = r⁵.
ϕphs(::Type{Val{5}}, x, xᵢ) = euclidean(x, xᵢ)^5
∂phs(::Type{Val{5}}, dim::Int) = ∂ℒ(x, xᵢ) = 5 * (x[dim] - xᵢ[dim]) * sqeuclidean(x, xᵢ)^1.5
∇phs(::Type{Val{5}}) = ∇ℒ(x, xᵢ) = 5 * (x .- xᵢ) * sqeuclidean(x, xᵢ)^1.5
function ∂²phs(::Type{Val{5}}, dim::Int)
    return function ∂²ℒ(x, xᵢ)
        return 5 * euclidean(x, xᵢ) * (3 * (x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ))
    end
end

# Septic polyharmonic spline, ϕ(r) = r⁷.
ϕphs(::Type{Val{7}}, x, xᵢ) = euclidean(x, xᵢ)^7
∂phs(::Type{Val{7}}, dim::Int) = ∂ℒ(x, xᵢ) = 7 * (x[dim] - xᵢ[dim]) * sqeuclidean(x, xᵢ)^2.5
∇phs(::Type{Val{7}}) = ∇ℒ(x, xᵢ) = 7 * (x .- xᵢ) * sqeuclidean(x, xᵢ)^2.5
function ∂²phs(::Type{Val{7}}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return 7 * sqeuclidean(x, xᵢ)^1.5 * (5 * (x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ))
    end
end

function Base.show(io::IO, rbf::PHS)
    return print(io, "Polyharmonic spline, rⁿ (n = $(rbf.n))")
end
