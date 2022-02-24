```@meta
CurrentModule = Pixell
```

# Working with Maps

## Plotting

We provide a simple plot recipe for displaying Enmap with Plots.jl. This package supports multiple backends: if you're using the default GR backend (i.e. by calling `using Plots; gr()`), then you can set the usual LaTeX font with

```julia
using Plots
default(fontfamily="Computer Modern", dpi=200)
```

The PyPlot backend uses matplotlib. If you're using this, i.e. using `using Plots; pyplot()`, then we recommend these settings,

```julia
using Plots
default(left_margin=15Plots.PlotMeasures.mm, dpi=200)
Plots.pyrcparams["text.usetex"] = true
Plots.pyrcparams["font.family"] = "serif"
Plots.pyrcparams["mathtext.fontset"] = "cm"
```

Finally, you can eschew our plot recipe completely and use `PyPlot.jl` to have a classic matplotlib experience. 

## Map Manipulation

## Relating pixels to the sky

