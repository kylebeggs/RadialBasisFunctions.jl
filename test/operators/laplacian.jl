using RadialBasisFunctions
using StaticArrays
using Statistics
using Random

rsme(test, correct) = sqrt(sum((test - correct) .^ 2) / sum(correct .^ 2))
mean_percent_error(test, correct) = mean(abs.((test .- correct) ./ correct)) * 100

f(x) = 1 + sin(4 * x[1]) + cos(3 * x[1]) + sin(2 * x[2])
d2f_dxx(x) = -16 * sin(4 * x[1]) - 9 * cos(3 * x[1])
d2f_dyy(x) = -4 * sin(2 * x[2])
∇²f(x) = d2f_dxx(x) + d2f_dyy(x)

N = 100
Δ = 1 / (N - 1)
points = 0:Δ:1
structured_points = ((x, y) for x in points for y in points)
x = map(x -> SVector{2}(x .+ (Δ / 5 .* rand(2))), structured_points)
y = f.(x)

@testset "Laplacian" begin
    ∇² = laplacian(x, PHS(3; poly_deg=4))
    @test mean_percent_error(∇²(y), ∇²f.(x)) < 5

    ∇² = laplacian(x, IMQ(1; poly_deg=4))
    @test mean_percent_error(∇²(y), ∇²f.(x)) < 5

    ∇² = laplacian(x, Gaussian(1; poly_deg=4))
    @test mean_percent_error(∇²(y), ∇²f.(x)) < 5
end

@testset "Different Evaluation Points" begin
    x2 = map(x -> SVector{2}(rand(MersenneTwister(x), 2)), (N + 1):(N + 1000))
    ∇² = laplacian(x, x2, PHS(3; poly_deg=4))
    @test mean_percent_error(∇²(y), ∇²f.(x2)) < 5
end
