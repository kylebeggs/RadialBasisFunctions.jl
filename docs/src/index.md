```@meta
CurrentModule = RadialBasisFunctions
```

# RadialBasisFunctions.jl

Documentation for [RadialBasisFunctions](https://github.com/kylebeggs/RadialBasisFunctions.jl).

This package intends to provide functions for all things regarding Radial Basis Functions (RBF). This includes

* Interpolation
* Linear Operators
  * Partial derivative ($\partial f$)
  * Laplacian ($\nabla^2 f$, $\Delta f$)
  * Gradient ($\nabla f$)
  * Directional Derivative ($\nabla f \cdot \mathbf{v}$)
  * Support for custom / user supplied operators ($\mathcal{L}$)

Currently, we support the following types of RBFs (all have polynomial augmentation by default, but is optional)

| Type                 | Function, $\phi(r)$                    |
| -------------------- | -------------------------------------- |
| Polyharmonic Spline  | $r^n$ where $n=1,3,5,7$          |
| Inverse Multiquadric | $1 / \sqrt{(r \varepsilon)^2+1}$ |
| Gaussian             | $e^{-(r \varepsilon)^2}$               |


where $\varepsilon$ is a user-supplied shape parameter.

## Installation

Simply install the latest stable release using Julia's package manager:

```julia
] add RadialBasisFunctions
```

## Planned Features

* Adaptive operators and interpolation. Adding / removing / modifying points and automatically updating the weights without a complete recalculation.
* Add more built-in operator combinations that will allow you to lazily construct operators such as
  * divergence
  * curl

## Current Limitations

1. A critical dependency of this package is [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) which requires that the dimension of each data point is inferrable. To quote from NearestNeighbors.jl:
    > The data, i.e., the points to build up the tree from. It can either be
    > * a matrix of size nd × np with the points to insert in the tree where nd is the dimensionality of the points and np is the number of points
    > * a vector of vectors with fixed dimensionality, nd, which must be part of the type. Specifically, data should be a Vector{V}, where V is itself a subtype of an AbstractVector and such that eltype(V) and length(V) are defined. (For example, with 3D points, V = SVector{3, Float64} works because eltype(V) = Float64 and length(V) = 3 are defined in V.)

    That said, we currently only support the second option here (`Vector{AbstractVector}`), but plan to support matrix inputs in the future.

2. `RadialBasisInterp` uses all points, but there are plans to support local collocation / subdomains like the operators use.
