"""
    abstract type AbstractRadialBasis end
"""
abstract type AbstractRadialBasis end

include("polyharmonic_spline.jl")
include("inverse_multiquadric.jl")
include("gaussian.jl")
include("monomials/monomial.jl")
include("monomials/monomial_builtin_operators.jl")

function ∂(basis::B, order::T, dim::T) where {T<:Int,B<:AbstractRadialBasis}
    return ∂(basis, Val{order}(), dim)
end
∂(basis::B, dim::T) where {T<:Int,B} = ∂(basis, Val{1}(), dim)
∂²(basis::B, dim::T) where {T<:Int,B} = ∂(basis, Val{2}(), dim)

# pretty printing
unicode_order(::Val{1}) = ""
unicode_order(::Val{2}) = "²"
unicode_order(::Val{3}) = "³"

function Base.show(io::IO, basis::B) where {B<:AbstractRadialBasis}
    if basis.poly_deg < 0
        print(io, "\n  No Monomial augmentation")
    else
        print(io, "\n  Monomial: degree $(basis.poly_deg)")
    end
end
