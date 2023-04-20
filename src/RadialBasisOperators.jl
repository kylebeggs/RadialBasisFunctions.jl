module RadialBasisOperators

using NearestNeighbors
using SymRCM
using LinearAlgebra
using LoopVectorization
using SparseArrays
using StaticArrays
using Statistics
using Distances
using ForwardDiff
using Enzyme
using StructArrays
using Enzyme
using Combinatorics
using FastClosures

abstract type AbstractRadialBasisOperator end

include("basis/basis.jl")
include("utils.jl")
include("linalg/stencil.jl")
include("operators/operators.jl")
include("interpolation.jl")

export RadialBasisOperator

# Scalar valued
export Partial, partial
export Laplacian, laplacian

# Vector valued
export Gradient, gradient

# interpolation
export RBFInterp

# Basis
export AbstractRadialBasis
export AbstractPHS, PHS, PHS1, PHS3, PHS5, PHS7
export IMQ
export Gaussian
export MonomialBasis
export ∂, ∂², ∇, ∇²

# utility functions
export find_neighbors, reorder_points!

# linear algebra
export _build_weightmx, _build_collocation_matrix!, _build_rhs!, _build_stencil!

#test
export ∂test, build_monomial_basis, pascals_triangle

end
