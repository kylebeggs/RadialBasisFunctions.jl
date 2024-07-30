using RadialBasisFunctions
using StaticArrays
using Statistics
using Random

rsme(test, correct) = sqrt(sum((test - correct) .^ 2) / sum(correct .^ 2))
mean_percent_error(test, correct) = mean(abs.((test .- correct) ./ correct)) * 100

f(x) = 1 + sin(4 * x[1]) + cos(3 * x[1]) + sin(2 * x[2])
df_dx(x) = 4 * cos(4 * x[1]) - 3 * sin(3 * x[1])
df_dy(x) = 2 * cos(2 * x[2])
d2f_dxx(x) = -16 * sin(4 * x[1]) - 9 * cos(3 * x[1])
d2f_dyy(x) = -4 * sin(2 * x[2])

N = 100
Δ = 1 / (N - 1)
points = 0:Δ:1
structured_points = ((x, y) for x in points for y in points)
x = map(x -> SVector{2}(x .+ (Δ / 2 .* rand(2))), structured_points)
y = f.(x)

@testset "First Derivative Partials" begin
    @testset "Polyharmonic Splines" begin
        ∂x = partial(x, 1, 1, PHS(3; poly_deg=2))
        ∂y = partial(x, 1, 2, PHS(3; poly_deg=2))
        @test mean_percent_error(∂x(y), df_dx.(x)) < 1
        @test mean_percent_error(∂y(y), df_dy.(x)) < 1
    end

    @testset "Inverse Multiquadrics" begin
        ∂x = partial(x, 1, 1, IMQ(1; poly_deg=2))
        ∂y = partial(x, 1, 2, IMQ(1; poly_deg=2))
        @test mean_percent_error(∂x(y), df_dx.(x)) < 1
        @test mean_percent_error(∂y(y), df_dy.(x)) < 1
    end

    @testset "Gaussian" begin
        ∂x = partial(x, 1, 1, Gaussian(1; poly_deg=2))
        ∂y = partial(x, 1, 2, Gaussian(1; poly_deg=2))
        @test mean_percent_error(∂x(y), df_dx.(x)) < 1
        @test mean_percent_error(∂y(y), df_dy.(x)) < 1
    end
end

@testset "Second Derivative Partials" begin
    @testset "Polyharmonic Splines" begin
        ∂2x = partial(x, 2, 1, PHS(3; poly_deg=4))
        ∂2y = partial(x, 2, 2, PHS(3; poly_deg=4))
        @test mean_percent_error(∂2x(y), d2f_dxx.(x)) < 1
        @test mean_percent_error(∂2y(y), d2f_dyy.(x)) < 1
    end

    @testset "Inverse Multiquadrics" begin
        ∂2x = partial(x, 2, 1, IMQ(1; poly_deg=4))
        ∂2y = partial(x, 2, 2, IMQ(1; poly_deg=4))
        @test mean_percent_error(∂2x(y), d2f_dxx.(x)) < 1
        @test mean_percent_error(∂2y(y), d2f_dyy.(x)) < 1
    end

    @testset "Gaussian" begin
        ∂2x = partial(x, 2, 1, Gaussian(1; poly_deg=4))
        ∂2y = partial(x, 2, 2, Gaussian(1; poly_deg=4))
        @test mean_percent_error(∂2x(y), d2f_dxx.(x)) < 1
        @test mean_percent_error(∂2y(y), d2f_dyy.(x)) < 1
    end
end

@testset "Different Evaluation Points" begin
    x2 = map(x -> SVector{2}(rand(2)), 1:100)
    ∂x = partial(x, x2, 1, 1, PHS(3; poly_deg=2))
    ∂y = partial(x, x2, 1, 2, PHS(3; poly_deg=2))
    @test mean_percent_error(∂x(y), df_dx.(x2)) < 1
    @test mean_percent_error(∂y(y), df_dy.(x2)) < 1
end
