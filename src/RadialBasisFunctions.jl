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

include("basis/basis.jl")
include("utils.jl")
include("linalg/stencil.jl")

include("operators/operators.jl")
include("operators/partial.jl")
include("operators/laplacian.jl")
include("operators/gradient.jl")
include("operators/monomial.jl")
include("operators/operator_combinations.jl")

include("interpolation.jl")

const Δ = ∇² # some people like this notation for the Laplacian
const DIV0_OFFSET = 1e-8

# operators
export RadialBasisOperator

# scalar valued
export Partial, partial
export Laplacian, laplacian

# vector valued
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
