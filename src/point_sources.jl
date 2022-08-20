using Pixell
using LinearAlgebra
using DelimitedFiles, Test
# import Pixell: plan_fftlog


##
using PythonCall
using Test
rft = RadialFourierTransform(n=256, pad=128)

pixell_utils = pyimport("pixell.utils")
np = pyimport("numpy")

ref_rft = pixell_utils.RadialFourierTransform(n=256, pad=128)
pyf = @pyeval `lambda x: 1/x`
h = real2harm(rft, r -> 1/r)
h_ref_1 = pyconvert(Vector, ref_rft.real2harm(pyf))
maximum(abs.(1 .- h ./ h_ref_1)) < 1000eps()


##
h = harm2real(rft, r -> 1/r)
h_ref_2 = pyconvert(Vector, ref_rft.harm2real(pyf))
# 1 .- h ./ h_ref
maximum(abs.(1 .- h ./ pyconvert(Vector, h_ref_2))) < 1000eps()

##
open("test/data/radialfouriertransform.txt", "w") do io
    writedlm(io, [h_ref_1 h_ref_2])
end
