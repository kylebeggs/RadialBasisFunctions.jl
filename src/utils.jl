function find_neighbors(data::AbstractVector, k::Int)
    tree = KDTree(data)
    adjl, _ = knn(tree, data, k, true)
    return adjl
end

function find_neighbors(data::AbstractVector, centers::AbstractVector, k::Int)
    tree = KDTree(data)
    adjl, _ = knn(tree, centers, k, true)
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
    x::AbstractVector{D}, adjl::Vector{Vector{T}}, k::T
) where {D,T<:Int}
    i = symrcm(adjl, ones(T, length(x)) .* k)
    permute!(x, i)
    return nothing
end

function reorder_points!(x::AbstractVector{D}, k::T) where {D,T<:Int}
    return reorder_points!(x, find_neighbors(x, k), k)
end

function check_poly_deg(poly_deg)
    if poly_deg < -1
        throw(ArgumentError("Augmented Monomial degree must be at least 0 (constant)."))
    end
    return nothing
end

_allocate_weights(m, n, k; sparse=true) = _allocate_weights(Float64, m, n, k; sparse=sparse)
function _allocate_weights(T, m, n, k; sparse=true)
    return sparse ? spzeros(T, m, n) : [zeros(T, k) for _ in 1:m]
end
