module RadialBasisOperators

using NearestNeighbors
using LinearAlgebra
using LoopVectorization
using SparseArrays
using StaticArrays
using Statistics
using Distances
using ForwardDiff
using Combinatorics
using Transducers

abstract type AbstractRadialBasisOperator end
abstract type VectorValuedOperator <: AbstractRadialBasisOperator end
abstract type ScalarValuedOperator <: AbstractRadialBasisOperator end

include("basis/basis.jl")
include("utils.jl")
include("linalg/linalg.jl")
include("operators/operators.jl")
include("interpolation.jl")

#export GenericRadialBasisOperator

# Scalar valued
export Partial, Partial², Laplacian

# Vector valued
export Gradient, Hessian

# Basis
export AbstractRadialBasis
export IMQ, PHS, Gaussian, polynomial_basis, ∂, ∂², ∇, ∇², find_neighbors
export PolynomialBasis

export build_weightmx, compute_stencil_weights
export RBFInterp, find_scales

end
