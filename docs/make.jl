using RadialBasisOperators
using Documenter

DocMeta.setdocmeta!(
    RadialBasisOperators, :DocTestSetup, :(using RadialBasisOperators); recursive=true
)

makedocs(;
    modules=[RadialBasisOperators],
    authors="Kyle Beggs",
    repo="https://github.com/kylebeggs/RadialBasisOperators.jl/blob/{commit}{path}#{line}",
    sitename="RadialBasisOperators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kylebeggs.github.io/RadialBasisOperators.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/kylebeggs/RadialBasisOperators.jl", devbranch="main")
