using RadialBasisFunctions
const RBF = RadialBasisFunctions
using StaticArrays
using LinearAlgebra

x = [SVector(1.0, 2.0), SVector(2.0, 1.0), SVector(1.5, 0.0)]

rb = PHS(3; poly_deg=1)
mb = MonomialBasis(2, 1)
L(x) = RBF.∂(x, 1, 1)
Lrb = L(rb)
Lmb = L(mb)

k = length(x)
n = k + 3

function unordered_approx(A::AbstractVector, B::AbstractVector)
    length(A) != length(B) && return false
    b = deepcopy(B)
    for a in A
        id = findfirst(x -> x ≈ a, b)
        isnothing(id) && return false
        deleteat!(b, id)
    end
    return true
end

@testset "Coefficient Matrix" begin
    A = Symmetric(zeros(n, n))
    _build_collocation_matrix!(A, x, rb, mb, k)
    @testset "RBFs" begin
        @test A[1, 2] ≈ (sqrt(sum((x[1] .- x[2]) .^ 2)))^3
        @test A[1, 3] ≈ (sqrt(sum((x[1] .- x[3]) .^ 2)))^3
        @test A[2, 3] ≈ (sqrt(sum((x[2] .- x[3]) .^ 2)))^3
    end
    @testset "Monomials" begin
        @test unordered_approx(A[1, 4:6], [1.0, x[1][1], x[1][2]])
        @test unordered_approx(A[2, 4:6], [1.0, x[2][1], x[2][2]])
        @test unordered_approx(A[3, 4:6], [1.0, x[3][1], x[3][2]])
    end
end

@testset "Right-hand side" begin
    b = zeros(n)
    _build_rhs!(b, Lrb, Lmb, x, rb, k)
    @testset "RBFs" begin
        @test b[1] ≈ Lrb(x[1], x[1])
        @test b[2] ≈ Lrb(x[1], x[2])
        @test b[3] ≈ Lrb(x[1], x[3])
    end
    @testset "Monomials" begin
        bb = zeros(3)
        Lmb(bb, x[1])
        @test all(b[4:6] .≈ bb)
    end
end
