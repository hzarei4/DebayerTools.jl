using DebayerTools
using Documenter

DocMeta.setdocmeta!(DebayerTools, :DocTestSetup, :(using DebayerTools); recursive=true)

makedocs(;
    modules=[DebayerTools],
    authors="Hossein Zarei Oshtolagh",
    sitename="DebayerTools.jl",
    format=Documenter.HTML(;
        canonical="https://hzarei4.github.io/DebayerTools.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/hzarei4/DebayerTools.jl",
    devbranch="master",
)
