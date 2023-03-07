########################################################################################
# Inverse Multiquadrics
struct IMQ{T,D<:Int} <: AbstractRadialBasis
    ε::T
    deg::D
    function IMQ(ε::T=1; deg::D=2) where {T,D<:Int}
        if all(ε .< 0)
            throw(ArgumentError("Shape parameter should be > 0. (ε=$ε)"))
        end
        deg < -1 && throw(
            ArgumentError("Augmented polynomial degree must be at least 0 (constant).")
        )
        return new{T,D}(ε, deg)
    end
end

(rbf::IMQ)(x, xᵢ) = 1 / sqrt((euclidean(x, xᵢ) * rbf.ε)^2 + 1)

function ∂(rbf::IMQ, dim::Int=1)
    function ∂ℒ(x, xᵢ)
        ε2 = rbf.ε .^ 2
        return (xᵢ[dim] - x[dim]) .* (ε2 / (ε2 * sqeuclidean(x, xᵢ) + 1)^1.5)
    end
end

function ∇(rbf::IMQ)
    function ∇ℒ(x, xᵢ)
        ε2 = rbf.ε .^ 2
        return (xᵢ - x) .* (ε2 / (ε2 * sqeuclidean(x, xᵢ) + 1)^1.5)
    end
end

function ∂²(rbf::IMQ, dim::Int=1)
    function ∂²ℒ(x, xᵢ)
        ε2 = rbf.ε .^ 2
        ε4 = ε2^2
        num1 = 3 * ε4 * (x[dim] - xᵢ[dim])^2
        denom = (ε2 * sqeuclidean(x, xᵢ) + 1)
        return num1 / denom^2.5 - ε2 / denom^1.5
    end
end

function franke(x)
    # modified Franke's formula for double precision as initial guess
    D = 2.0 * euclidean(first(x), last(x))
    N = size(points, 2)
    return D / (0.8 * (N^0.25))
end

function Base.show(io::IO, rbf::IMQ)
    return print(io, "Inverse Multiquadrics, 1/sqrt((r*ε)²+1) (ε = $(rbf.ε))")
end
