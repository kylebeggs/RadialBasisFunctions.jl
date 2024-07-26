using RadialBasisFunctions
import RadialBasisFunctions as RBF
using StaticArrays

@testset "dim=1, deg=0, singleton input" begin
    x = 2

    m = MonomialBasis(1, 0)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{1,0}

    # standard evaluation
    @test m(x) ≈ 1

    # derivatives
    @test RBF.∂(m, 1)(x) ≈ 0
    @test RBF.∂²(m, 1)(x) ≈ 0
    @test RBF.∇²(m)(x) ≈ 0
    grad = RBF.∇(m)(x)
    @test grad[1] ≈ 0
end

@testset "dim=1, deg=0, vector input" begin
    x = SVector(2.0)

    m = MonomialBasis(1, 0)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{1,0}

    # standard evaluation
    @test all(isapprox.(m(x), [1]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0]))
    @test all(isapprox.(RBF.∇²(m)(x), [0]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0]))
end

@testset "dim=1, deg=1" begin
    x = SVector(2.0)

    m = MonomialBasis(1, 1)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{1,1}

    # standard evaluation
    @test all(isapprox.(m(x), [1, 2]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0, 1]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0, 0]))
    @test all(isapprox.(RBF.∇²(m)(x), [0, 0]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0, 1]))
end

@testset "dim=1, deg=2" begin
    x = SVector(2.0)

    m = MonomialBasis(1, 2)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{1,2}

    # standard evaluation
    @test all(isapprox.(m(x), [1, 2, 4]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0, 1, 4]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0, 0, 2]))
    @test all(isapprox.(RBF.∇²(m)(x), [0, 0, 2]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0, 1, 4]))
end

@testset "dim=2, deg=0" begin
    x = SVector(2.0, 3.0)

    m = MonomialBasis(2, 0)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{2,0}

    # standard evaluation
    @test all(isapprox.(m(x), [1]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0]))
    @test all(isapprox.(RBF.∂(m, 2)(x), [0]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0]))
    @test all(isapprox.(RBF.∂²(m, 2)(x), [0]))
    @test all(isapprox.(RBF.∇²(m)(x), [0]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0]))
    @test all(isapprox.(grad[2], [0]))
end

@testset "dim=2, deg=1" begin
    x = SVector(2.0, 3.0)

    m = MonomialBasis(2, 1)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{2,1}

    # standard evaluation
    @test all(isapprox.(m(x), [1, 2, 3]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0, 1, 0]))
    @test all(isapprox.(RBF.∂(m, 2)(x), [0, 0, 1]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0, 0, 0]))
    @test all(isapprox.(RBF.∂²(m, 2)(x), [0, 0, 0]))
    @test all(isapprox.(RBF.∇²(m)(x), [0, 0, 0]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0, 1, 0]))
    @test all(isapprox.(grad[2], [0, 0, 1]))
end

@testset "dim=2, deg=2" begin
    x = SVector(2.0, 3.0)

    m = MonomialBasis(2, 2)
    @test m isa MonomialBasis
    @test typeof(m) <: MonomialBasis{2,2}

    # standard evaluation
    @test all(isapprox.(m(x), [1, 2, 3, 6, 4, 9]))

    # derivatives
    @test all(isapprox.(RBF.∂(m, 1)(x), [0, 1, 0, 3, 4, 0]))
    @test all(isapprox.(RBF.∂(m, 2)(x), [0, 0, 1, 2, 0, 6]))
    @test all(isapprox.(RBF.∂²(m, 1)(x), [0, 0, 0, 0, 2, 0]))
    @test all(isapprox.(RBF.∂²(m, 2)(x), [0, 0, 0, 0, 0, 2]))
    @test all(isapprox.(RBF.∇²(m)(x), [0, 0, 0, 0, 2, 2]))
    grad = RBF.∇(m)(x)
    @test all(isapprox.(grad[1], [0, 1, 0, 3, 4, 0]))
    @test all(isapprox.(grad[2], [0, 0, 1, 2, 0, 6]))
end
