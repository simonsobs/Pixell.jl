module Pixell

using WCS
using FITSIO
using FFTW
using Printf

include("enmap.jl")
include("enmap_ops.jl")

export Enmap, CarClenshawCurtis
export fullsky_geometry
export read_map, write_map

end
