"""
    PolynomialBasis{N<:Int,D<:Int,B}

Multivariate polynomial basis. Adapted from https://discourse.julialang.org/t/generation-of-multivariate-polynomial-basis-with-given-degree/60757/3?
n ∈ N: length of array, i.e., x ∈ Rⁿ
d ∈ N: degree
"""
struct PolynomialBasis{N<:Int,D<:Int}
    n::N
    d::D
    basis::Function
    function PolynomialBasis(n::N, deg::D) where {N<:Int,D<:Int}
        exponents = multiexponents(n + 1, deg)
        function basis(x::T) where {T}
            x = vcat(x, one(eltype(x)))
            return collect(Map(exponent -> prod(x .^ exponent))(exponents))
        end
        return new{N,D}(n, deg, basis)
    end
end

(p::PolynomialBasis)(x) = p.basis(x)

∂(p::PolynomialBasis, dim::Int) = ∂ℒ(x) = ForwardDiff.jacobian(p.basis, x)[:, dim]
∇(p::PolynomialBasis) = ∇ℒ(x) = ForwardDiff.jacobian(p.basis, x)
∂²(p::PolynomialBasis, dim::Int) = ∂²ℒ(x) = ForwardDiff.jacobian(∂(p.basis, dim), x)[:, dim]
∇²(p::PolynomialBasis) = ∇ℒ(x) = sum(ForwardDiff.hessian(p.basis, x), 2)
