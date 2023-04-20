using RadialBasisOperators
using StaticArrays
using LinearAlgebra

x = [SVector(1.0, 2.0), SVector(2.0, 1.0), SVector(1.5, 0.0)]

rb = PHS(3; poly_deg=1)
mb = MonomialBasis(2, 1)
L(x) = ∂(x, Val{1}(), 1)
Lrb = L(rb)
Lmb = L(mb)

k = length(x)
n = k + 3

A = Symmetric(zeros(n, n))
b = zeros(n)

_build_collocation_matrix!(A, x, rb, mb, k)
_build_rhs!(b, Lrb, Lmb, x, rb, k)

@testset "Coefficient Matrix" begin
    @testset "RBFs" begin
        @test A[1, 2] ≈ (sqrt(sum((x[1] .- x[2]) .^ 2)))^3
        @test A[1, 3] ≈ (sqrt(sum((x[1] .- x[3]) .^ 2)))^3
        @test A[2, 3] ≈ (sqrt(sum((x[2] .- x[3]) .^ 2)))^3
    end
    @testset "Monomials" begin
        @test all(A[1, 4:6] .≈ (1.0, x[1][1], x[1][2]))
        @test all(A[2, 4:6] .≈ (1.0, x[2][1], x[2][2]))
        @test all(A[3, 4:6] .≈ (1.0, x[3][1], x[3][2]))
    end
end

@testset "Right-hand side" begin
    @testset "RBFs" begin
        @test b[1] ≈ Lrb(x[1], x[1])
        @test b[2] ≈ Lrb(x[1], x[2])
        @test b[3] ≈ Lrb(x[1], x[3])
    end
    @testset "Monomials" begin
        @test all(b[4:6] .≈ Lmb(x[1]))
    end
end
