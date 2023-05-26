# Getting Started

First, let's load the package along with the [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package which we use for each data point

```@example overview
using RadialBasisFunctions
using StaticArrays
```

## Interpolation

Suppose we have a set of data ``\mathbf{x}`` where ``\mathbf{x}_i \in \mathbb{R}^2``, and we want to interpolate a function ``f:\mathbb{R}^2 \rightarrow \mathbb{R}``

```@example overview
f(x) = 2*x[1]^2 + 3*x[2]
x = [SVector{2}(rand(2)) for _ in 1:1000]
y = f.(x)
```

and now we can build the interpolator

```@example overview
interp = RadialBasisInterp(x, y)
```

and evaluate it at a new point

```@example overview
x_new = [rand(2) for _ in 1:5]
y_new = interp(x_new)
y_true = f.(x_new)
```

and compare the error

```@example overview
abs.(y_true .- y_new)
```

Wow! The error is numerically zero! Well... we set ourselves up for success here. `RadialBasisInterp` (along with `RadialBasisOperator`) has an optional argument to provide the type of radial basis including the degree of polynomial augmentation. The default basis is a cubic polyharmonic spline with 2nd degree polynomial augmentation (which the constructor is `PHS(3, poly_deg=2)`) and given the underlying function we are interpolating is a 2nd order polynomial itself, we are able to represent it exactly (up to machine precision). Let's see what happens when we only use 1st order polynomial augmentation

```@example overview
interp = RadialBasisInterp(x, y, PHS(3, poly_deg=1))
y_new = interp(x_new)
abs.(y_true .- y_new)
```

## Operators

more to come here...
