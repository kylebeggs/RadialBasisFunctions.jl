########################################################################################
# Polyharmonic Spline

"""
   abstract type AbstractPHS <: AbstractRadialBasis

Supertype of all Polyharmonic Splines.
"""
abstract type AbstractPHS <: AbstractRadialBasis end

# convienience constructor
function PHS(n::T=3; poly_deg::T=2) where {T<:Int}
    check_poly_deg(poly_deg)
    if iseven(n) || n > 7
        throw(ArgumentError("n must be 1, 3, 5, or 7. (n = $n)"))
    end
    n == 1 && return PHS1(poly_deg)
    n == 3 && return PHS3(poly_deg)
    n == 5 && return PHS5(poly_deg)
    return PHS7(poly_deg)
end

# Linear polyharmonic spline, ϕ(r) = r.
struct PHS1{T<:Int} <: AbstractPHS
    poly_deg::T
    function PHS1(poly_deg::T) where {T<:Int}
        check_poly_deg(poly_deg)
        return new{T}(poly_deg)
    end
end

(phs::PHS1)(x, xᵢ) = euclidean(x, xᵢ)
∂(::PHS1, ::Val{1}, dim::Int) = ∂ℒ(x, xᵢ) = (x[dim] - xᵢ[dim]) / euclidean(x, xᵢ)
∇(::PHS1) = ∇ℒ(x, xᵢ) = (x .- xᵢ) / euclidean(x, xᵢ)
function ∂(::PHS1, ::Val{2}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return (-(x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ)) / (euclidean(x, xᵢ)^3 + 1e-8)
    end
end
function ∇²(::PHS1)
    function ∇²ℒ(x, xᵢ)
        return sum((-(x .- xᵢ) .^ 2 .+ sqeuclidean(x, xᵢ)) / (euclidean(x, xᵢ)^3 + 1e-8))
    end
end

# Cubic polyharmonic spline, ϕ(r) = r³.
struct PHS3{T<:Int} <: AbstractPHS
    poly_deg::T
    function PHS3(poly_deg::T) where {T<:Int}
        check_poly_deg(poly_deg)
        return new{T}(poly_deg)
    end
end

(phs::PHS3)(x, xᵢ) = euclidean(x, xᵢ)^3
∂(::PHS3, ::Val{1}, dim::Int) = ∂ℒ(x, xᵢ) = 3 * (x[dim] - xᵢ[dim]) * euclidean(x, xᵢ)
∇(::PHS3) = ∇ℒ(x, xᵢ) = 3 * (x .- xᵢ) * euclidean(x, xᵢ)
function ∂(::PHS3, ::Val{2}, dim::Int)
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
struct PHS5{T<:Int} <: AbstractPHS
    poly_deg::T
    function PHS5(poly_deg::T) where {T<:Int}
        check_poly_deg(poly_deg)
        return new{T}(poly_deg)
    end
end

(phs::PHS5)(x, xᵢ) = euclidean(x, xᵢ)^5
∂(::PHS5, ::Val{1}, dim::Int) = ∂ℒ(x, xᵢ) = 5 * (x[dim] - xᵢ[dim]) * euclidean(x, xᵢ)^3
∇(::PHS5) = ∇ℒ(x, xᵢ) = 5 * (x .- xᵢ) * euclidean(x, xᵢ)^3
function ∂(::PHS5, ::Val{2}, dim::Int)
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
struct PHS7{T<:Int} <: AbstractPHS
    poly_deg::T
    function PHS7(poly_deg::T) where {T<:Int}
        check_poly_deg(poly_deg)
        return new{T}(poly_deg)
    end
end

(phs::PHS7)(x, xᵢ) = euclidean(x, xᵢ)^7
∂(::PHS7, ::Val{1}, dim::Int) = ∂ℒ(x, xᵢ) = 7 * (x[dim] - xᵢ[dim]) * euclidean(x, xᵢ)^5
∇(::PHS7) = ∇ℒ(x, xᵢ) = 7 * (x .- xᵢ) * euclidean(x, xᵢ)^5
function ∂(::PHS7, ::Val{2}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        return 7 * euclidean(x, xᵢ)^3 * (5 * (x[dim] - xᵢ[dim])^2 + sqeuclidean(x, xᵢ))
    end
end
function ∇²(::PHS7)
    function ∇²ℒ(x, xᵢ)
        return sum(7 * euclidean(x, xᵢ)^3 * (5 * (x .- xᵢ) .^ 2 .+ sqeuclidean(x, xᵢ)))
    end
end

function Base.show(io::IO, rbf::R) where {R<:AbstractPHS}
    print(io, print_basis(rbf))
    if rbf.poly_deg < 0
        print(io, "\n  No Monomial augmentation")
    else
        print(io, "\n└─Polynomial augmentation: degree $(rbf.poly_deg)")
    end
end

print_basis(::PHS1) = "Polyharmonic spline (r¹)"
print_basis(::PHS3) = "Polyharmonic spline (r³)"
print_basis(::PHS5) = "Polyharmonic spline (r⁵)"
print_basis(::PHS7) = "Polyharmonic spline (r⁷)"
