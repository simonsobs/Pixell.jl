```@meta
CurrentModule = Pixell
```

# Developer Notes

If you're new to the [Julia language](https://julialang.org/), my favorite tutorial is this [Introduction to Julia for Quantitative Economics](https://julia.quantecon.org/intro.html). Sections 1-13 are the most relevant if you are an astronomer, but some of later numerical techniques can also be useful. Computing clusters will often have a Julia module; make sure you load a Julia version 1.6 or later. 

For your development computer, download the [precompiled binaries](https://julialang.org/downloads/) and add the Julia executable to your `$PATH`. We **strongly** advise against using a package manager to install Julia. The Julia compiler has a heavily modified glibc (similar to Rust), which many package managers misconfigure (such as the AUR). 


## Package development

To make changes to this package, start up the Julia interpeter and run

```julia-repl
julia> ] dev git@github.com:simonsobs/Pixell.jl.git
```

By default, this will place a copy of this package's repository in your home directory, `~/.julia/dev/Pixell`. Changes to the code in this folder will be reflected in your global environment.

It can be helpful to write documentation with live preview. This is accomplished with the [LiveServer](https://github.com/tlienart/LiveServer.jl) package. If you have a standard Julia installlation, you can run from the command line,

```bash
julia --project=$HOME/.julia/dev/Pixell/docs -e "using Pixell, LiveServer; servedocs()"
```
This will render HTML pages and provide an HTTP server for local previews. The documentation will automatically update when you make changes.
