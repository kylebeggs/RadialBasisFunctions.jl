function find_neighbors(data::AbstractVector, k::Int)
    tree = KDTree(data)
    adjl, _ = knn(tree, data, k, true)
    return adjl
end

function find_neighbors(data::AbstractVector, eval_points::AbstractVector, k::Int)
    tree = KDTree(data)
    adjl, _ = knn(tree, eval_points, k, true)
    return adjl
end

function get_num_params(f)
    return length.(getfield.(getfield.(methods(f), :sig), :parameters)) .- 1
end

function check_num_params(ℒ, basis::B) where {B<:AbstractRadialBasis}
    num_params_ϕ = get_num_params(ℒ(basis.ϕ))
    !any(num_params_ϕ .== 2) && throw(ArgumentError("ϕ must have 2 input arguments."))
    #num_params_poly = get_num_params(ℒ(basis.poly))
    #!any(num_params_poly .== 1) && throw(ArgumentError("poly must have 1 input argument."))
    return nothing
end

function check_if_deg_odd(deg::Int)
    if isodd(deg)
        @warn "Monomial degree is recommended to be even, unless it is -1 which indicates no Monomials. (Flyer, 2016 - https://doi.org/10.1016/j.jcp.2016.05.026)"
    end
end

"""
    autoselect_k(data::Vector, basis<:AbstractRadialBasis)

See Bayona, 2017 - https://doi.org/10.1016/j.jcp.2016.12.008
"""
function autoselect_k(data::Vector, basis::B) where {B<:AbstractRadialBasis}
    m = basis.poly_deg
    d = length(first(data))
    return min(length(data), max(2 * binomial(m + d, d), 2 * d + 1))
end

function reorder_points!(
    x::AbstractVector, adjl::AbstractVector{AbstractVector{T}}, k::T
) where {T<:Int}
    i = symrcm(adjl, ones(T, length(x)) .* k)
    permute!(x, i)
    return nothing
end

function reorder_points!(x::AbstractVector, k::T) where {T<:Int}
    return reorder_points!(x, find_neighbors(x, k), k)
end

function check_poly_deg(poly_deg)
    if poly_deg < -1
        throw(ArgumentError("Augmented Monomial degree must be at least 0 (constant)."))
    end
    return nothing
end

function scale_cloud(data)
    furthest_point = maximum(p -> euclidean(first(data), p), data)
    return data ./ furthest_point
end

_allocate_weights(m, n, k) = _allocate_weights(Float64, m, n, k)
function _allocate_weights(T, m, n, k)
    return spzeros(T, m, n)
end

function columnwise_div(A::SparseMatrixCSC, B::AbstractVector)
    I, J, V = findnz(A)
    for idx in eachindex(V)
        V[idx] /= B[I[idx]]
    end
    return sparse(I, J, V)
end
columnwise_div(A::SparseMatrixCSC, B::Number) = A ./ B

function _find_smallest_dist(data, k)
    tree = KDTree(data)
    _, dists = knn(tree, data, k, true)
    Δ = minimum(dists) do d
        z = minimum(@view(d[2:end])) do i
            abs(i - first(d))
        end
        return z
    end
    return Δ
end

_get_underlying_type(x::AbstractVector) = eltype(x)
_get_underlying_type(x::Number) = typeof(x)
