```@meta
CurrentModule = Pixell
```

# Pixell.jl

[Pixell](https://github.com/simonsobs/Pixell.jl) is a package for high performance data analysis of sky maps with rectangular pixels. It is based on a subset of the Python package [pixell](https://github.com/simonsobs/pixell). This package has a particular focus on astrophysical and cosmological science enabled by efficient map manipulation, easy multithreading and GPU support, and machine learning. 

This package manipulates maps on equidistant cylindrical projections (ECP) like [plate carrée](https://pro.arcgis.com/en/pro-app/2.8/help/mapping/properties/plate-carree.htm). It implements an array type equipped with WCS information. Another common pixelization of the sphere is [Healpix.jl](https://github.com/ziotom78/Healpix.jl), which uses a constant size pixel with a more complicated pixel shape. 

Pixell supports development and deployment on a wide variety of [platforms](https://julialang.org/downloads/), ranging from laptops to computing clusters. Installation is as simple as starting up the Julia interpreter, and running

```julia-repl
julia> ] add Pixell
```
