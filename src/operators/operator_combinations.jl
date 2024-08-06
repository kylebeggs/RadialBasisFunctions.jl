for op in (:+, :-, :*, :/)
    @eval function Base.$op(a::ℒRadialBasisFunction, b::ℒRadialBasisFunction)
        return additive_ℒRBF(x, xᵢ) = Base.$op(a(x, xᵢ), b(x, xᵢ))
    end
end

for op in (:+, :-, :*, :/)
    @eval function Base.$op(a::ℒMonomialBasis, b::ℒMonomialBasis)
        function additive_ℒMon(m, x)
            m .= Base.$op.(a(x), b(x))
            return nothing
        end
    end
end

for op in (:+, :-, :*, :/)
    @eval function Base.$op(op1::RadialBasisOperator, op2::RadialBasisOperator)
        if !all(op1.data .≈ op2.data)
            throw(
                ArgumentError(
                    "Can not add operators that were not built with the same data."
                ),
            )
        end
        if !all(op1.adjl .≈ op2.adjl)
            throw(
                ArgumentError("Can not add operators that do not have the same stencils.")
            )
        end
        k1 = length(first((op1.adjl)))
        k2 = length(first((op2.adjl)))
        k = k1 > k2 ? k1 : k2
        ℒ(x) = Base.$op(op1.ℒ(x), op2.ℒ(x))
        return RadialBasisOperator(ℒ, op1.data, op1.basis; k=k, adjl=op1.adjl)
    end
end
