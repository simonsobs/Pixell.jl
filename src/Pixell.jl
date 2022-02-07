module Pixell

using WCS
using FITSIO
using FFTW
using Printf
import Unitful, UnitfulAngles
import Unitful: uconvert, ustrip

include("enmap.jl")
include("enmap_geom.jl")
include("enmap_ops.jl")

export Enmap, CarClenshawCurtis, getwcs
export geometry, fullsky_geometry, slice_geometry
export pix2sky, pix2sky!, sky2pix, sky2pix!
export read_map, write_map

# set up some shortcuts for common angles
const radian = Unitful.rad
const degree = Unitful.°
const arcminute = UnitfulAngles.arcminute

end
