struct RadialBasisInterp{X,Y,R,M,RB,MB}
    x::X
    y::Y
    rbf_weights::R
    monomial_weights::M
    rbf_basis::RB
    monomial_basis::MB
end

function RadialBasisInterp(x, y, basis::B=PHS()) where {B<:AbstractRadialBasis}
    dim = length(first(x))
    k = length(x)  # number of data in influence/support domain
    npoly = binomial(dim + basis.poly_deg, basis.poly_deg)
    n = k + npoly
    poly = MonomialBasis(dim, basis.poly_deg)
    A = Symmetric(zeros(eltype(first(x)), n, n))
    _build_collocation_matrix!(A, x, basis, k, poly)
    b = zeros(eltype(first(x)), n)
    b[1:k] .= y
    w = A \ b
    return RadialBasisInterp(x, y, w[1:k], w[(k + 1):end], basis, poly)
end

function (rbfi::RadialBasisInterp)(x::T) where {T}
    rbf = zero(eltype(T))
    for i in eachindex(rbfi.rbf_weights)
        rbf += rbfi.rbf_weights[i] * rbfi.rbf_basis(x, rbfi.x[i])
    end

    val_poly = rbfi.monomial_basis(x)
    poly = zero(eltype(T))
    for i in eachindex(rbfi.monomial_weights)
        poly += rbfi.monomial_weights[i] * val_poly[i]
    end
    return rbf + poly
end

(rbfi::RadialBasisInterp)(x::Vector{<:AbstractVector}) = [rbfi(val) for val in x]
