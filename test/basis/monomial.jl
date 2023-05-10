using RadialBasisFunctions
const RBF = RadialBasisFunctions
using StaticArrays

function unordered_approx(A::AbstractVector, B::AbstractVector)
    length(A) != length(B) && return false
    b = deepcopy(B)
    for a in A
        id = findfirst(x -> x ≈ a, b)
        isnothing(id) && return false
        deleteat!(b, id)
    end
    return true
end

x = SVector(2.0, 3.0)

m = MonomialBasis(2, 2)
@test m isa MonomialBasis
@test m.n == 2
@test m.deg == 2

# standard evaluation
@test unordered_approx(m(x), [1, 2, 3, 4, 9, 6])

# derivatives
dx! = RBF.∂(m, 1, 1)
dy! = RBF.∂(m, 1, 2)
d2x! = RBF.∂(m, 2, 1)
d2y! = RBF.∂(m, 2, 2)
grad = RBF.∇(m)
lap! = RBF.∇²(m)

b = zeros(6)

dx!(b, x)
@test unordered_approx(b, [0, 1, 0, 4, 0, 3])

dy!(b, x)
@test unordered_approx(b, [0, 0, 1, 0, 6, 2])

d2x!(b, x)
@test unordered_approx(b, [0, 0, 0, 2, 0, 0])

d2y!(b, x)
@test unordered_approx(b, [0, 0, 0, 0, 2, 0])

lap!(b, x)
@test unordered_approx(b, [0, 0, 0, 2, 2, 0])

b = ntuple(x -> zeros(6), length(x))
grad(b, x)
@test unordered_approx(b[1], [0, 1, 0, 4, 0, 3])
@test unordered_approx(b[2], [0, 0, 1, 0, 6, 2])
