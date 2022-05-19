
using PyCall
enmap, curvedsky = pyimport("pixell.enmap"), pyimport("pixell.curvedsky")

##
using PowerSpectra
using Pixell
shape, wcs = fullsky_geometry(20 * Pixell.arcminute)

Dl = PowerSpectra.planck_theory_Dl()
channelindex(X) = findfirst(first(X), "TEB")

ells = collect(eachindex(Dl[:TT]))
Dlfac = ells .* (ells .+1 ) ./ (2π)

𝐂 = zeros(3, 3, length(ells))
for X in ("T", "E", "B"), Y in ("T", "E", "B")
    c₁, c₂ = channelindex(X), channelindex(Y)
    XY = Symbol(X * Y)
    if XY in keys(Dl)
        𝐂[c₁, c₂ ,:] .= parent(Dl[XY]) ./ Dlfac
        𝐂[ c₂, c₁ ,:] .= 𝐂[c₁, c₂ ,:]
    end
end

lmax = Pixell.getlmax(wcs)
alms = (Alm(lmax, lmax), Alm(lmax, lmax), Alm(lmax, lmax))
synalm!(𝐂, alms; lmin=2, lmax=lmax)

##
m = alm2map(alms[1], shape, wcs)
plot(m)

##
I, Q, U = alm2map(alms, shape, wcs)
aT, aE, aB = map2alm((I, Q, U); lmax=lmax)
ClTT = alm2cl(aT)
ClEE = alm2cl(aE)
ClTE = alm2cl(aT, aE)
ClBB = alm2cl(aB)

IQU = cat(I, Q, U, dims=3)
write_map("data/IQU.fits", IQU)

##
aT.lmax
##
using DelimitedFiles

open("data/TEB_alms_real.dat", "w") do io
    writedlm(io, [real.(a.alm) for a in alms])
end;
open("data/TEB_alms_imag.dat", "w") do io
    writedlm(io, [imag.(a.alm) for a in alms])
end;


##
pyI = pycall(enmap.read_map, PyObject, "data/IQU.fits", sel=PyObject(0:0))  # get WCS
newpymap = pycall(enmap.enmap, PyObject, pyI, pyI.wcs)
pyI = pycall(curvedsky.alm2map, PyObject, PyObject(alms[1].alm), newpymap, spin=0)

pyQU = pycall(enmap.read_map, PyObject, "data/IQU.fits", sel=PyObject(1:2))
newpymap = pycall(enmap.enmap, PyObject, pyQU, pyQU.wcs)
pyQUalm = [alms[2].alm  alms[3].alm]
pyQU = pycall(curvedsky.alm2map, PyObject, PyReverseDims(pyQUalm), newpymap, spin=2)


enmap.write_map("data/pyI.fits", pyI)
enmap.write_map("data/pyQU.fits", pyQU)

pyI = pycall(enmap.read_map, PyObject, "data/IQU.fits", sel=PyObject(0:0))
pyQU = pycall(enmap.read_map, PyObject, "data/IQU.fits", sel=PyObject(1:2))

# islice = py"lambda x : x[3,:,:]"
# quslice = py"lambda x : x[3,:,:]"

# pyI = pycall(islice, PyObject, pyIQU)

##
pyT = curvedsky.map2alm(pyI, lmax=lmax, spin=0)
pyEB = curvedsky.map2alm(pyQU, lmax=lmax, spin=2)
pyE = pyEB[1,:]
pyB = pyEB[2,:]

aT, aE, aB = map2alm((I, Q,U); lmax=lmax)
hp = pyimport("healpy")

pyClTT = (hp.alm2cl(pyT)')[:,1]
pyClTE = (hp.alm2cl(pyT, pyE)')[:,1]
pyClEE = hp.alm2cl(pyE)
pyClBB = hp.alm2cl(pyB)

##
using DelimitedFiles

data = [pyClTT pyClTE pyClEE pyClBB]
open("data/test_cls_IQU.txt", "w") do io
    writedlm(io, data)
end;
