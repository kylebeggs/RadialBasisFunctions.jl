function find_neighbors(data::Vector, k::Int)
    tree = KDTree(data)
    adjl, _ = knn(tree, data, k, true)
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
        @warn "Polynomial degree is recommended to be even, unless it is -1 which indicates no polynomials. (Flyer, 2016 - https://doi.org/10.1016/j.jcp.2016.05.026)"
    end
end

"""
    autoselect_k(data::Vector, basis<:AbstractRadialBasis)

See Bayona, 2017 - https://doi.org/10.1016/j.jcp.2016.12.008
"""
function autoselect_k(data::Vector, basis::B) where {B<:AbstractRadialBasis}
    m = basis.deg
    d = length(first(data))
    return min(length(data), max(2 * binomial(m + d, d), 2 * d + 1))
end

function find_scales(data::Vector, adjl::Vector{<:AbstractVector})
    scales = map(adjl) do adj
        stencil = @view data[adj]
        cur_max = typemin(eltype(eltype((stencil))))
        for i in eachindex(stencil), j in eachindex(stencil)
            dist = euclidean(stencil[i], stencil[j])
            dist > cur_max && (cur_max = dist)
        end
        return cur_max
    end
    return scales
end
