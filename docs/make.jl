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
        assets=["assets/custom.css"]
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Map Manipulation" => "map_manipulation.md",
        "Developer Notes" => "developer_notes.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/simonsobs/Pixell.jl",
    devbranch="main",
)
