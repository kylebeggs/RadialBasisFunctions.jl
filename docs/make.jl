using RadialBasisFunctions
using Documenter

DocMeta.setdocmeta!(
    RadialBasisFunctions, :DocTestSetup, :(using RadialBasisFunctions); recursive=true
)

makedocs(;
    modules=[RadialBasisFunctions],
    authors="Kyle Beggs",
    repo="https://github.com/kylebeggs/RadialBasisFunctions.jl/blob/{commit}{path}#{line}",
    sitename="RadialBasisFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kylebeggs.github.io/RadialBasisFunctions.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/kylebeggs/RadialBasisFunctions.jl", devbranch="main")
