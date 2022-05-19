
# this is a work-in-progress test involving a picture of a duck

using Pixell, Plots
duck = read_map("data/duck.fits")
shape, wcs = size(duck), duck.wcs
plot(duck)

##
alm_duck = map2alm(duck)
nm = Pixell.alm2map(alm_duck, size(duck), duck.wcs)
plot(nm)

##
right_duck = duck[end:-1:begin, :]

plot(right_duck)

duckQU = cat(duck, right_duck; dims=3)
write_map("data/duckQU.fits", duckQU)

##
almDuck_QU = map2alm(QU)
newduck = alm2map(almDuck_QU, shape, wcs)
plot(newduck[1] .- duck)

##

pyQU = pycall(enmap.read_map, PyObject, "data/duckQU.fits")

myalm = pycall(curvedsky.map2alm, PyObject, pyQU, lmax=almDuck_QU[1].lmax, spin=2)
pyE = collect(myalm)[1,:][1]
pyB = collect(myalm)[2,:][1]

myalm2 = [pyE pyB]

##
almDuck_QU[1].alm
plot(abs.(almDuck_QU[1].alm), yscale=:log10, ylim=(1e-10,1e-1), lw=0.1)
plot!(abs.(pyE), yscale=:log10, ylim=(1e-10,1e-1), lw=0.1)

##
newpymap = pycall(enmap.enmap, PyObject, pyQU, pyQU.wcs)
# pyQUalm = collect(cat(almDuck_QU[1].alm,almDuck_QU[2].alm,dims=2)')
pyQUduck = pycall(curvedsky.alm2map, PyObject, PyReverseDims(myalm2), newpymap, spin=2)

##
l2 = py"lambda x : x[1,:,:]"
heatmap(l2(pyQUduck))

##
# heatmap(duck.data .- nm.data, clim=(-1e-5,1e-5))
#
# imap = pycall(enmap.read_map, PyObject, "test/data/duck.fits")
# pyshape = imap.shape
# pywcs = imap.wcs
# pya = curvedsky.map2alm(imap, lmax=alm_duck.lmax)
# newpymap = pycall(enmap.enmap, PyObject, imap,  imap.wcs)
# pycall(curvedsky.alm2map, PyObject, pya, newpymap, spin=0)
# heatmap(PyArray(imap)' - PyArray(newpymap)', clim=(-1e-5,1e-5))
