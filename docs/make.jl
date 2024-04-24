using RadialBasisFunctions
using Documenter

DocMeta.setdocmeta!(
    RadialBasisFunctions, :DocTestSetup, :(using RadialBasisFunctions); recursive=true
)

makedocs(;
    modules=[RadialBasisFunctions],
    authors="Kyle Beggs",
    sitename="RadialBasisFunctions.jl",
    repo = Documenter.Remotes.GitHub("kylebeggs", "RadialBasisFunctions.jl"),
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kylebeggs.github.io/RadialBasisFunctions.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
)

withenv("GITHUB_REPOSITORY" => repo) do
  deploydocs(; repo, versions=["stable" => "v^", "dev" => "dev"])
end
