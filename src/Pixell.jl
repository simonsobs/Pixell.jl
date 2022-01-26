module Pixell

using WCS
using FITSIO
using FFTW
using Printf

include("enmap.jl")
include("enmap_ops.jl")

export Enmap, CarClenshawCurtis, getwcs
export fullsky_geometry, slice_geometry
export pix2sky, pix2sky!, sky2pix, sky2pix!
export read_map, write_map

end
