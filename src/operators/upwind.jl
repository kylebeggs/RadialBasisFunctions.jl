struct Upwind{L<:Function,T<:Int} <: AbstractOperator
    ℒ::L
    dim::T
end

# convienience constructors
"""
    function upwind(data, dim, basis; k=autoselect_k(data, basis))

Builds a `RadialBasisOperator` where the operator is the partial derivative, `Partial`, of `order` with respect to `dim`.
"""
function upwind(
    data::AbstractVector{D},
    dim,
    basis::B=PHS(3; poly_deg=2);
    k::T=autoselect_k(data, basis),
) where {D<:AbstractArray,T<:Int,B<:AbstractRadialBasis}
    tree = KDTree(data)
    _, dists = knn(tree, data, k, true)
    Δ = minimum(dists) do d
        z = minimum(@view(d[2:end])) do i
            abs(i - first(d))
        end
        return z
    end

    N = length(first(data))
    dx = zeros(N)
    dx[dim] = Δ

    li = regrid(data, data .- Ref(dx); k=k)
    ri = regrid(data, data .+ Ref(dx); k=k)
    update_weights!(li)
    update_weights!(ri)
    one_typed = one(eltype(first(data)))
    l = columnwise_div((sparse(one_typed * I, size(li.weights)...) .- li.weights), Δ)
    r = columnwise_div((sparse(one_typed * I, size(ri.weights)...) .- ri.weights), Δ)
    c = partial(data, 1, dim; k=k)

    du = let l = l, r = r, c = c
        (ϕ, v, θ) -> begin
            wl = max.(v, Ref(0)) .* (θ * (l * ϕ) .+ (1 - θ) * c(ϕ))
            wr = min.(v, Ref(0)) .* (θ * (r * ϕ) .+ (1 - θ) * c(ϕ))
            wl .+ wr
        end
    end

    return du
end

function columnwise_div(A::SparseMatrixCSC, B::AbstractVector)
    I, J, V = findnz(A)
    for idx in eachindex(V)
        V[idx] /= B[I[idx]]
    end
    return sparse(I, J, V)
end
columnwise_div(A::SparseMatrixCSC, B::Number) = A ./ B

separate(f, x) = (t = findall(f, x); (@view(x[t]), @view(x[setdiff(eachindex(x), t)])))

function split_stencils(x, adjl, dim)
    split = map(adjl) do a
        center = first(a)
        l, r = separate(i -> x[i][dim] < x[center][dim], @view(a[2:end]))
        return vcat(center, l), vcat(center, r)
    end
    return getindex.(split, 1), getindex.(split, 2)
end
