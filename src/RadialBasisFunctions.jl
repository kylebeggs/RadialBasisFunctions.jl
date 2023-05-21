module RadialBasisFunctions

using NearestNeighbors
using SymRCM
using LinearAlgebra
using LoopVectorization
using SparseArrays
using StaticArrays
using Statistics
using Distances
using Combinatorics

abstract type AbstractRadialBasisOperator end

include("basis/basis.jl")
include("utils.jl")
include("linalg/stencil.jl")
include("operators/operators.jl")
include("interpolation.jl")

const Δ = ∇² # some people like this notation for the Laplacian

# operators
export RadialBasisOperator

# Scalar valued
export Partial, partial
export Laplacian, laplacian

# Vector valued
export Gradient, gradient

# interpolation
export RadialBasisInterp

# bases
export AbstractRadialBasis
export RadialBasisFunction
export AbstractPHS, PHS, PHS1, PHS3, PHS5, PHS7
export IMQ
export Gaussian
export MonomialBasis

# utility functions
export find_neighbors, reorder_points!

# linear algebra
export _build_weightmx, _build_collocation_matrix!, _build_rhs!, _build_stencil!

#test
export ∂test, ∂exponents, build_monomial_basis, pascals_triangle, monomial_recursive_list

end