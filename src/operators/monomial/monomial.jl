include("partial.jl")
#include("gradient.jl")
#include("laplacian.jl")

function ∂(::MonomialBasis{Dim,Deg,O}, order::T, dim::T) where {Dim,Deg,T<:Int,O}
    if order == 1
        return MonomialBasis{Dim,Deg}(∂_get_ops(Val(Dim), Val(Deg), Val(dim)))
    elseif order == 2
        return MonomialBasis{Dim,Deg}(∂_get_ops(Val(Dim), Val(Deg), Val(dim)))
    else
        throw(
            ArgumentError("Only first and second order partial derivatives are supported")
        )
    end
end

∂²(mb::MonomialBasis, dim::T) where {T} = ∂(mb, 2, dim)

function ∇(m::MonomialBasis)
    ∂! = ntuple(dim -> ∂(m, 1, dim), m.n)
    function ∇ℒ(b, x)
        for i in eachindex(∂!)
            b[i] .= 0
            ∂![i](b[i], x)
        end
        return nothing
    end
    return ∇ℒ
end

function ∇²(m::MonomialBasis)
    ∂² = ntuple(dim -> ∂(m, 2, dim), m.n)
    function ∇²ℒ(b, x)
        cache = ones(size(b))
        b .= 0
        for ∂²! in ∂²
            # use mapreduce here instead?
            ∂²!(cache, x)
            b .+= cache
        end
        return nothing
    end
    return ∇²ℒ
end
