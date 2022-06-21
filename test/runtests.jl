using Pixell
using Test
import Pixell: degree, arcminute

using DelimitedFiles

include("test_geometry.jl")  # creating geometries and sky ↔ pix
include("test_enmap.jl")     # enmap features and manipulation
include("test_transforms.jl")
include("test_distance_transform.jl")
include("test_io.jl")
include("test_plot.jl")
include("test_utils.jl")