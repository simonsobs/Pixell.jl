```@meta
CurrentModule = Pixell
```

# Tutorial

Let's make an Enmap, the primary structure in this package.

```@example tutorial
using Plots # hide
Plots.default(fontfamily="Computer Modern", fmt=:svg) # hide
```

```@example tutorial
using Pixell, Plots
Plots.default(fontfamily="Computer Modern", fmt=:svg)
shape, wcs = fullsky_geometry(300.0 * Pixell.arcminute)  # set up the map geometry
m = Enmap(randn(shape), wcs)  # generate a random map with the shape and WCS
plot(m)
```

Let's compute a spherical harmonic transform with Libsharp, and then compute the power spectrum ``C_{\ell}``.
```@example tutorial
cl = alm2cl(map2alm(m))
plot(cl, ylabel=raw"$C_{\ell}$", xlabel=raw"Multipole moment, $\ell$")
```


## Reading and writing maps

## Making some plots


## Noteworthy differences from Python
The Julia language has a number of [differences from Python](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python), and the [Pixell.jl](https://github.com/simonsobs/Pixell.jl) has a number of important differences from the Python package [pixell](https://github.com/simonsobs/pixell).
