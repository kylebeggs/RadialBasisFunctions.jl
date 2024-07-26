"""
    abstract type AbstractBasis end
"""
abstract type AbstractBasis end

"""
    abstract type AbstractRadialBasis <: AbstractBasis end
"""
abstract type AbstractRadialBasis <: AbstractBasis end

struct ℒRadialBasisFunction{F<:Function}
    f::F
end
(ℒrbf::ℒRadialBasisFunction)(x, xᵢ) = ℒrbf.f(x, xᵢ)

struct ℒMonomial{F<:Function}
    f::F
end
(ℒmon::ℒMonomial)(x) = ℒmon.f(x)

include("polyharmonic_spline.jl")
include("inverse_multiquadric.jl")
include("gaussian.jl")
include("monomial.jl")

# pretty printing
unicode_order(::Val{1}) = ""
unicode_order(::Val{2}) = "²"
unicode_order(::Val{3}) = "³"
unicode_order(::Val{4}) = "⁴"
unicode_order(::Val{5}) = "⁵"
unicode_order(::Val{6}) = "⁶"
unicode_order(::Val{7}) = "⁷"
unicode_order(::Val{8}) = "⁸"
unicode_order(::Val{9}) = "⁹"

function Base.show(io::IO, basis::B) where {B<:AbstractRadialBasis}
    if basis.poly_deg < 0
        print(io, "\n  No Monomial augmentation")
    else
        print(io, "\n  Monomial: degree $(basis.poly_deg)")
    end
end
