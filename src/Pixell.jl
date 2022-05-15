module Pixell

using WCS
using WCS_jll
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
import Healpix: Alm, map2alm, alm2cl

include("enmap.jl")
include("enmap_geom.jl")
include("enmap_ops.jl")
include("transforms.jl")
include("plot.jl")
include("transform_distance.jl")

export Enmap, CarClenshawCurtis, getwcs
export geometry, fullsky_geometry, slice_geometry
export pix2sky, pix2sky!, sky2pix, sky2pix!, skyarea, pixareamap, pixareamap!
export read_map, write_map
export Alm, map2alm, alm2cl
export distance_transform, ApproxSeqSDT, ExactSeqSDT, BruteForceSDT

# set up some shortcuts for common angles
const radian = Unitful.rad
const degree = Unitful.°
const arcminute = UnitfulAngles.arcminute

function __init__()
    Enplot.register_colors!()
end


end
