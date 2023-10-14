module RadialBasisFunctions

using NearestNeighbors
using SymRCM
using LinearAlgebra
using LoopVectorization
using ChunkSplitters
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

using PrecompileTools
@setup_workload begin
    f(x) = 1 + sin(4 * x[1]) + cos(3 * x[1]) + sin(2 * x[2])
    x = map(x -> SVector{2}(rand(2)), 1:100)
    z = f.(x)
    @compile_workload begin
        # basis functions
        basis_funcs = [
            IMQ(1),
            Gaussian(1),
            PHS(1; poly_deg=-1),
            PHS(3; poly_deg=0),
            PHS(5; poly_deg=1),
            PHS(7; poly_deg=2),
        ]

        for b in basis_funcs
            # partial
            ∂ = partial(x, 1, 1, b)
            ∂x = ∂(z)

            # gradient
            ∇ = gradient(x, b)
            ∇z = ∇(z)

            # laplacian
            ∇² = laplacian(x, b)
            ∇²z = ∇²(z)

            # interpolation
            interp = RadialBasisInterp(x, z, b)
            zz = interp([SVector{2}(rand(2)), SVector{2}(rand(2))])
        end
    end
end

end # module
