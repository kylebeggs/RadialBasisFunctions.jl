########################################################################################
# Polyharmonic Spline

abstract type AbstractPHS <: AbstractRadialBasis end

"""
    PHS

Polyharmonic spline , ϕ(r) = rⁿ. Allowable values of n are 1, 3, 5, or 7.
"""
function PHS(n::T=3; poly_deg::T=2) where {T<:Int}
    iseven(n) && throw(ArgumentError("n should be odd (n=$n)."))
    poly_deg < -1 &&
        throw(ArgumentError("Augmented Monomial degree must be at least 0 (constant)."))
    return eval(Symbol("PHS" * string(n)))(poly_deg)
end

# Linear polyharmonic spline, ϕ(r) = r.
struct PHS1{D<:Int} <: AbstractPHS
    poly_deg::D
end

(phs::PHS1)(x, xᵢ) = euclidean(x, xᵢ)
∂(::PHS1, dim::Int) = ∂ℒ(x, xᵢ) = (x[dim] - xᵢ[dim]) / euclidean(x, xᵢ)
∇(::PHS1) = ∇ℒ(x, xᵢ) = (x .- xᵢ) / euclidean(x, xᵢ)
function ∂²(::PHS1, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return (-(x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ)) /
               (sqeuclidean(x, xᵢ)^1.5 + 1e-8)
    end
end
function ∇²(::PHS1)
    function ∇²ℒ(x, xᵢ)
        return sum(
            (-(x .- xᵢ) .^ 2 .+ sqeuclidean(x, xᵢ)) / (sqeuclidean(x, xᵢ)^1.5 + 1e-8)
        )
    end
end

# Cubic polyharmonic spline, ϕ(r) = r³.
struct PHS3{D<:Int} <: AbstractPHS
    poly_deg::D
end

(phs::PHS3)(x, xᵢ) = euclidean(x, xᵢ)^3
∂(::PHS3, dim::Int) = ∂ℒ(x, xᵢ) = 3 * (x[dim] - xᵢ[dim]) * euclidean(x, xᵢ)
∇(::PHS3) = ∇ℒ(x, xᵢ) = 3 * (x .- xᵢ) * euclidean(x, xᵢ)
function ∂²(::PHS3, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return 3 * (sqeuclidean(x, xᵢ) + (x[dim] - xᵢ[dim])^2) / (euclidean(x, xᵢ) + 1e-8)
    end
end
function ∇²(::PHS3)
    function ∇²ℒ(x, xᵢ)
        return sum(3 * (sqeuclidean(x, xᵢ) .+ (x .- xᵢ) .^ 2) / (euclidean(x, xᵢ) + 1e-8))
    end
end

# Quintic polyharmonic spline, ϕ(r) = r⁵.
struct PHS5{D<:Int} <: AbstractPHS
    poly_deg::D
end

(phs::PHS5)(x, xᵢ) = euclidean(x, xᵢ)^5
∂(::PHS5, dim::Int) = ∂ℒ(x, xᵢ) = 5 * (x[dim] - xᵢ[dim]) * sqeuclidean(x, xᵢ)^1.5
∇(::PHS5) = ∇ℒ(x, xᵢ) = 5 * (x .- xᵢ) * sqeuclidean(x, xᵢ)^1.5
function ∂²(::PHS5, dim::Int)
    return function ∂²ℒ(x, xᵢ)
        return 5 * euclidean(x, xᵢ) * (3 * (x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ))
    end
end
function ∇²(::PHS5)
    return function ∇²ℒ(x, xᵢ)
        return sum(5 * euclidean(x, xᵢ) * (3 * (x .- xᵢ) .^ 2 .+ sqeuclidean(x, xᵢ)))
    end
end

# Septic polyharmonic spline, ϕ(r) = r⁷.
struct PHS7{D<:Int} <: AbstractPHS
    poly_deg::D
end

(phs::PHS7)(x, xᵢ) = euclidean(x, xᵢ)^7
∂(::PHS7, dim::Int) = ∂ℒ(x, xᵢ) = 7 * (x[dim] - xᵢ[dim]) * sqeuclidean(x, xᵢ)^2.5
∇(::PHS7) = ∇ℒ(x, xᵢ) = 7 * (x .- xᵢ) * sqeuclidean(x, xᵢ)^2.5
function ∂²(::PHS7, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return 7 * sqeuclidean(x, xᵢ)^1.5 * (5 * (x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ))
    end
end
function ∇²(::PHS7)
    function ∇²ℒ(x, xᵢ)
        return sum(7 * sqeuclidean(x, xᵢ)^1.5 * (5 * (x .- xᵢ) .^ 2 .+ sqeuclidean(x, xᵢ)))
    end
end

function Base.show(io::IO, rbf::R) where {R<:AbstractPHS}
    print(io, print_phs(rbf))
    if rbf.poly_deg < 0
        print(io, "\n  No Monomial augmentation")
    else
        print(io, "\n└─Polynomial augmentation: degree $(rbf.poly_deg)")
    end
end

print_phs(::PHS1) = "Polyharmonic spline, r¹"
print_phs(::PHS3) = "Polyharmonic spline, r³"
print_phs(::PHS5) = "Polyharmonic spline, r⁵"
print_phs(::PHS7) = "Polyharmonic spline, r⁷"
