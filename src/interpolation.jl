struct RBFInterp{X,Y,T,P}
    x::Vector{X}
    y::Vector{Y}
    rbf_weights::Vector{T}
    poly_weights::Vector{T}
    basis::AbstractRadialBasis
    poly::P
end

function RBFInterp(x, y, basis::B=PHS(3, 1)) where {B<:AbstractRadialBasis}
    dim = length(first(x))
    np = length(x)  # number of data in influence/support domain
    npoly = binomial(dim + basis.deg, basis.deg)
    poly = polynomial_basis(dim, basis.deg)
    C = build_collocation_matrix(x, basis, np, npoly, poly)
    b = zeros(np + npoly)
    b[1:np] = y
    w = C \ b
    return RBFInterp(x, y, w[1:np], w[(np + 1):end], basis, poly)
end

function (rbfi::RBFInterp)(x::T) where {T}
    rbf = zero(eltype(T))
    @tturbo for i in eachindex(rbfi.rbf_weights)
        rbf += rbfi.rbf_weights[i] * rbfi.basis(x, rbfi.x[i])
    end

    val_poly = rbfi.poly(x)
    poly = zero(eltype(T))
    @tturbo for i in eachindex(rbfi.poly_weights)
        poly += rbfi.poly_weights[i] * val_poly[i]
    end
    return rbf + poly
end

(rbfi::RBFInterp)(x::Vector{<:AbstractVector}) = [rbfi(val) for val in x]
