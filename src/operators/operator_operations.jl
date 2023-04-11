function Base.:+(a::Partial, b::Partial)
    if !all(a.data .≈ b.data)
        throw(
            ArgumentError("Can not add operators that were not built with the same data.")
        )
    end
    if !all(a.adjl .≈ b.adjl)
        throw(ArgumentError("Can not add operators that do not have the same stencils."))
    end
    order = vcat(a.order, b.order)
    dim = vcat(a.dim, b.dim)
    weights = a.weights .+ b.weights
    return partial(order, dim, weights, a.data, a.adjl, a.basis)
end
