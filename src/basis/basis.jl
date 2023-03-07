"""
    abstract type AbstractRadialBasis end
"""
abstract type AbstractRadialBasis end

include("polyharmonic_spline.jl")
include("inverse_multiquadric.jl")
include("gaussian.jl")
include("polynomial.jl")

function Base.show(io::IO, basis::B) where {B<:AbstractRadialBasis}
    print(io, "AbstractRadialBasis")
    print(io, "\n  Radial: ", basis)
    if basis.deg < 0
        print(io, "\n  No polynomial augmentation")
    else
        print(io, "\n  Polynomial: degree $(basis.deg)")
    end
end
