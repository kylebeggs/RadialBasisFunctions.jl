"""
    struct Gaussian{T,D<:Int} <: AbstractRadialBasis

Gaussian radial basis function:``ϕ(r) = e^{-(ε r)^2}``
"""
struct Gaussian{T,D<:Int} <: AbstractRadialBasis
    ε::T
    poly_deg::D
    function Gaussian(ε::T=1; poly_deg::D=2) where {T,D<:Int}
        if all(ε .< 0)
            throw(ArgumentError("Shape parameter should be > 0. (ε=$ε)"))
        end
        poly_deg < -1 &&
            throw(ArgumentError("Augmented Monomial degree must be at least 0 (constant)."))
        return new{T,D}(ε, poly_deg)
    end
end

(rbf::Gaussian)(x, xᵢ) = exp(-(rbf.ε * euclidean(x, xᵢ))^2)

function ∂(rbf::Gaussian, ::Val{1}, dim::Int)
    function ∂ℒ(x, xᵢ)
        return -2 * rbf.ε^2 * (x[dim] - xᵢ[dim]) * exp(-rbf.ε^2 * sqeuclidean(x, xᵢ))
    end
end

function ∇(rbf::Gaussian)
    function ∇ℒ(x, xᵢ)
        return -2 * rbf.ε^2 * (x .- xᵢ) * exp(-rbf.ε^2 * sqeuclidean(x, xᵢ))
    end
end

function ∂(rbf::Gaussian, ::Val{2}, dim::Int)
    function ∂²ℒ(x, xᵢ)
        ε2 = rbf.ε^2
        return (4 * ε2^2 * (x[dim] - xᵢ[dim])^2 - 2 * ε2) * exp(-ε2 * sqeuclidean(x, xᵢ))
    end
end

function ∇²(rbf::Gaussian)
    function ∇²ℒ(x, xᵢ)
        ε2 = rbf.ε^2
        return sum((4 * ε2^2 * (x .- xᵢ) .^ 2 .- 2 * ε2) * exp(-ε2 * sqeuclidean(x, xᵢ)))
    end
end

function Base.show(io::IO, rbf::Gaussian)
    print(io, "Gaussian, exp(-(ε*r)²)")
    print(io, "\n  ├─Shape factor: ε = $(rbf.ε)")
    if rbf.poly_deg < 0
        print(io, "\n  └─No Monomial augmentation")
    else
        print(io, "\n  └─Polynomial augmentation: degree $(rbf.poly_deg)")
    end
end

print_basis(rbf::Gaussian) = "Gaussian (ε = $(rbf.ε))"
