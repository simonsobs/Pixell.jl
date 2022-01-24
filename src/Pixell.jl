module Pixell

using WCS
using FITSIO
using FFTW
using Printf

import WCS: pix_to_world  # we extend this to Enmap

include("enmap.jl")
include("enmap_ops.jl")

export Enmap, CarClenshawCurtis
export fullsky_geometry
export pix_to_world

end
