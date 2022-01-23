using Pixell
using Test

using PyCall
enmap = pyimport("pixell.enmap")

@testset "Pixell.jl" begin
    # Write your tests here.
    shape, wcs_py = enmap.geometry(shape=(128, 256), res=deg2rad(5/60.),pos=(0,0))
    m = pycall(enmap.ones, PyObject, shape, wcs_py)
end
