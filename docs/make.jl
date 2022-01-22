using Pixell
using Documenter

DocMeta.setdocmeta!(Pixell, :DocTestSetup, :(using Pixell); recursive=true)

makedocs(;
    modules=[Pixell],
    authors="Zack Li, Yilun Guan",
    repo="https://github.com/xzackli/Pixell.jl/blob/{commit}{path}#{line}",
    sitename="Pixell.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/Pixell.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xzackli/Pixell.jl",
    devbranch="main",
)
