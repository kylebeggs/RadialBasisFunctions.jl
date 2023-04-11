using RadialBasisOperators
using StaticArrays
using Statistics

rsme(test, correct) = sqrt(sum((test - correct) .^ 2) / sum(correct .^ 2))
percent_error(test, correct) = mean(abs.((test .- correct) ./ correct)) * 100

f(x) = 1 + sin(4 * x[1]) + cos(3 * x[1]) + sin(2 * x[2])
df_dx(x) = 4 * cos(4 * x[1]) - 3 * sin(3 * x[1])
df_dy(x) = 2 * cos(2 * x[2])
d2f_dxx(x) = -16 * sin(4 * x[1]) - 9 * cos(3 * x[1])
d2f_dyy(x) = -4 * sin(2 * x[2])

N = 1000
x = map(x -> SVector{2}(rand(2)), 1:N)
y = f.(x)

@testset "First Derivative Partials" begin
    ∂x = partial(x, 1, 1, PHS(3; poly_deg=4))
    ∂y = partial(x, 1, 2, PHS(3; poly_deg=4))
    @test percent_error(∂x(y), df_dx.(x)) < 2
    @test percent_error(∂y(y), df_dy.(x)) < 2

    ∂x = partial(x, 1, 1, IMQ(1; poly_deg=4))
    ∂y = partial(x, 1, 2, IMQ(1; poly_deg=4))
    @test percent_error(∂x(y), df_dx.(x)) < 2
    @test percent_error(∂y(y), df_dy.(x)) < 2

    ∂x = partial(x, 1, 1, Gaussian(1; poly_deg=4))
    ∂y = partial(x, 1, 2, Gaussian(1; poly_deg=4))
    @test percent_error(∂x(y), df_dx.(x)) < 2
    @test percent_error(∂y(y), df_dy.(x)) < 2
end

@testset "Second Derivative Partials" begin
    ∂2x = partial(x, 2, 1, PHS(3; poly_deg=4))
    ∂2y = partial(x, 2, 2, PHS(3; poly_deg=4))
    @test percent_error(∂2x(y), d2f_dxx.(x)) < 2
    @test percent_error(∂2y(y), d2f_dyy.(x)) < 2

    ∂2x = partial(x, 2, 1, IMQ(1; poly_deg=4))
    ∂2y = partial(x, 2, 2, IMQ(1; poly_deg=4))
    @test percent_error(∂2x(y), d2f_dxx.(x)) < 2
    @test percent_error(∂2y(y), d2f_dyy.(x)) < 2

    ∂2x = partial(x, 2, 1, Gaussian(1; poly_deg=4))
    ∂2y = partial(x, 2, 2, Gaussian(1; poly_deg=4))
    @test percent_error(∂2x(y), d2f_dxx.(x)) < 2
    @test percent_error(∂2y(y), d2f_dyy.(x)) < 2
end
