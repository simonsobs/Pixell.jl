module Pixell

using WCS
import WCS: AbstractWCSTransform
using FITSIO
using FFTW
using Printf
import Unitful, UnitfulAngles
import Unitful: uconvert, ustrip
using StaticArrays
using DSP: unwrap, unwrap!
import FastTransforms: chebyshevjacobimoments1, clenshawcurtisweights
using Libsharp
import Libsharp: sharp_execute!
import Healpix: Alm, map2alm

include("enmap.jl")
include("enmap_geom.jl")
include("enmap_ops.jl")
include("transforms.jl")

export Enmap, CarClenshawCurtis, getwcs
export geometry, fullsky_geometry, slice_geometry
export pix2sky, pix2sky!, sky2pix, sky2pix!, skyarea
export read_map, write_map
export Alm, map2alm

# set up some shortcuts for common angles
const radian = Unitful.rad
const degree = Unitful.°
const arcminute = UnitfulAngles.arcminute

include("plot.jl")

end
