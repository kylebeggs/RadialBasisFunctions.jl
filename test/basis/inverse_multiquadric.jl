using RadialBasisFunctions
using StaticArrays

@testset "Constructors" begin
    imq = IMQ()
    @test imq isa IMQ
    @test imq.ε == 1
    @test imq.poly_deg == 2

    imq = IMQ(5.0; poly_deg=0)
    @test imq.ε ≈ 5
    @test imq.poly_deg == 0

    @test_throws ArgumentError IMQ(-1)
    @test_throws ArgumentError IMQ(; poly_deg=-2)
end

x₁ = SVector(1.0, 2)
x₂ = SVector(2.0, 4)
imq = IMQ(2; poly_deg=-1)

@testset "Distances" begin
    r = sqrt((x₁[1] - x₂[1])^2 + (x₁[2] - x₂[2])^2)
    @test imq(x₁, x₂) ≈ 1 / sqrt((imq.ε * r)^2 + 1)
end

@testset "Derivatives" begin
    dim = 1
    ∂rbf = ∂(imq, dim)
    ∂²rbf = ∂²(imq, dim)
    ∇rbf = ∇(imq)

    @test ∂rbf(x₁, x₂) ≈ 4 / (21 * sqrt(21))
    @test all(∇rbf(x₁, x₂) .≈ (4 / (21 * sqrt(21)), 8 / (21 * sqrt(21))))
    @test ∂²rbf(x₁, x₂) ≈ -4 / (49 * sqrt(21))
end
