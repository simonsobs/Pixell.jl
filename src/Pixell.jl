module Pixell

import Base: in
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
import Healpix: Alm, map2alm, alm2map, alm2cl
import AbstractFFTs: fftfreq, Plan, plan_fft!, plan_ifft!
import LinearAlgebra: mul!, ldiv!
import SpecialFunctions: loggamma

include("enmap.jl")

include("projections/arbitrary_wcs.jl")
include("projections/car_proj.jl")
include("projections/tan_proj.jl")

include("enmap_geom.jl")
include("enmap_ops.jl")
include("transforms.jl")
include("plot.jl")
include("transform_distance.jl")
include("utils.jl")

export Enmap, CarClenshawCurtis, Gnomonic, getwcs
export geometry, fullsky_geometry, slice_geometry, pad, SkyBoundingBox
export pix2sky, pix2sky!, sky2pix, sky2pix!, skyarea, pixareamap, pixareamap!
export read_map, write_map
export Alm, map2alm, alm2cl, alm2map, alm2map!
export distance_transform, ApproxSeqSDT, ExactSeqSDT, BruteForceSDT
export dplanck, RadialFourierTransform, real2harm, harm2real
export FFTLogPlan

# set up some shortcuts for common angles
const radian = Unitful.rad
const degree = Unitful.°
const arcminute = UnitfulAngles.arcminute

function __init__()
    Enplot.register_colors!()
end

end
