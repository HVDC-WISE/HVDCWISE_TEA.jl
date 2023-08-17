using HVDCWISE_TEA
using Documenter

DocMeta.setdocmeta!(HVDCWISE_TEA, :DocTestSetup, :(using HVDCWISE_TEA); recursive=true)

makedocs(;
    modules=[HVDCWISE_TEA],
    authors="Matteo Baù, Matteo Rossini",
    repo="https://github.com/HVDC-WISE/HVDCWISE_TEA.jl/blob/{commit}{path}#{line}",
    sitename="HVDCWISE_TEA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HVDC-WISE.github.io/HVDCWISE_TEA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HVDC-WISE/HVDCWISE_TEA.jl",
    devbranch="main",
)
