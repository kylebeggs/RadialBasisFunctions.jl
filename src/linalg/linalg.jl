
function build_weightmx(ℒ, data, adjl, basis, scales)
    dim = length(first(data))
    np = length(first(adjl))  # number of data in influence/support domain
    npoly = binomial(dim + basis.deg, basis.deg)
    poly = PolynomialBasis(dim, basis.deg)
    ℒpoly = ℒ(poly)
    num_neighbors = length.(adjl)
    I = zeros(Int, sum(num_neighbors))
    J = reduce(vcat, adjl)
    V = zeros(eltype(first(data)), sum(num_neighbors))
    Threads.@threads for i in eachindex(num_neighbors)
        range = ((i - 1) * np + 1):(i * np)
        @turbo I[range] .= i
        #d = data[adjl[i]] ./ scales[i]
        d = data[adjl[i]]
        V[range] = @views compute_stencil_weights(ℒ, d, basis, np, npoly, poly, ℒpoly)
    end
    return sparse(I, J, V)
end

function compute_stencil_weights(
    ℒ, data, basis::B, num_points, num_poly, poly, ℒpoly
) where {B<:AbstractRadialBasis}
    C = build_collocation_matrix(data, basis, num_points, num_poly, poly)
    b = build_rhs(ℒ, data, basis, num_points, num_poly, ℒpoly)
    return (C \ b)[1:num_points]
end

function build_collocation_matrix(data, basis, num_points, num_poly, poly)
    # radial basis section
    C = zeros(num_points + num_poly, num_points + num_poly)
    @inbounds for j in 1:num_points, i in 1:j
        C[i, j] = basis(data[i], data[j])
    end
    # add polynomial augmentation
    if basis.deg > -1
        @inbounds for i in 1:num_points
            C[i, (num_points + 1):end] = poly(data[i])
        end
    end
    return Symmetric(C, :U)
end

function build_rhs(ℒ, data, basis, num_points, num_poly, ℒpoly)
    b = zeros(num_points + num_poly)
    @inbounds for i in eachindex(data)
        b[i] = ℒ(basis)(first(data), data[i])
    end
    if basis.deg > -1
        b[(num_points + 1):end] = ℒpoly(first(data))
    end
    return b
end
