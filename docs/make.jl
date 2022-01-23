using Pixell
using Documenter

DocMeta.setdocmeta!(Pixell, :DocTestSetup, :(using Pixell); recursive=true)

makedocs(;
    modules=[Pixell],
    authors="Zack Li, Yilun Guan",
    repo="https://github.com/simonsobs/Pixell.jl/blob/{commit}{path}#{line}",
    sitename="Pixell.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://simonsobs.github.io/Pixell.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/simonsobs/Pixell.jl",
    devbranch="main",
)
