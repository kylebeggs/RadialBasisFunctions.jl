"""
    abstract type AbstractRadialBasis end
"""
abstract type AbstractRadialBasis end

include("polyharmonic_spline.jl")
include("inverse_multiquadric.jl")
include("gaussian.jl")
include("monomial.jl")

∂(basis::B, ::Type{Val{1}}, dim::T) where {T<:Int,B} = ∂(basis, dim)
∂(basis::B, ::Type{Val{2}}, dim::T) where {T<:Int,B} = ∂²(basis, dim)

unicode_order(::Type{Val{1}}) = ""
unicode_order(::Type{Val{2}}) = "²"
unicode_order(::Type{Val{3}}) = "³"

function Base.show(io::IO, basis::B) where {B<:AbstractRadialBasis}
    if basis.poly_deg < 0
        print(io, "\n  No Monomial augmentation")
    else
        print(io, "\n  Monomial: degree $(basis.poly_deg)")
    end
end
