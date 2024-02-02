using RadialBasisFunctions
using StaticArrays
using Statistics
using Random
using LinearAlgebra

rsme(test, correct) = sqrt(sum((test - correct) .^ 2) / sum(correct .^ 2))
mean_percent_error(test, correct) = mean(abs.((test .- correct) ./ correct)) * 100

f(x) = 1 + sin(4 * x[1]) + cos(3 * x[1]) + sin(2 * x[2])
df_dx(x) = 4 * cos(4 * x[1]) - 3 * sin(3 * x[1])
df_dy(x) = 2 * cos(2 * x[2])
d2f_dxx(x) = -16 * sin(4 * x[1]) - 9 * cos(3 * x[1])
d2f_dyy(x) = -4 * sin(2 * x[2])

N = 10_000
x = map(x -> SVector{2}(rand(MersenneTwister(x), 2)), 1:N)
y = f.(x)

@testset "Single Direction" begin
    v = SVector(2.0, 1.0)
    v /= norm(v)
    ∇v = directional(x, v, PHS3(2))
    ∇vy = ∇v(y)
    @test mean_percent_error(∇vy[1], df_dx.(x) .* v[1]) < 1
    @test mean_percent_error(∇vy[2], df_dy.(x) .* v[2]) < 1
end

@testset "Direction Vector for Each Data Center" begin
    v = map(1:length(x)) do i
        v = SVector{2}(rand(2))
        return v /= norm(v)
    end
    ∇v = directional(x, v, PHS3(2))
    ∇vy = ∇v(y)
    df_dx_v = map((df, v) -> df * v[1], df_dx.(x), v)
    df_dy_v = map((df, v) -> df * v[2], df_dy.(x), v)
    @test mean_percent_error(∇vy[1], df_dx_v) < 1
    @test mean_percent_error(∇vy[2], df_dy_v) < 1
end
