using SafeTestsets

@safetestset "Polyharmonic Splines" begin
    include("basis/polyharmonic_spline.jl")
end

@safetestset "Gaussian" begin
    include("basis/gaussian.jl")
end

@safetestset "Inverse Multiquadric" begin
    include("basis/inverse_multiquadric.jl")
end
