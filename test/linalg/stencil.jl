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
∇²f(x) = func_dxx(x) + func_dyy(x)

N = 1000
x = map(x -> SVector{2}(rand(2)), 1:N)
y = f.(x)
