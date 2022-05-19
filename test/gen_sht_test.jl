
using Pixell, Plots
duck = read_map("test/data/duck.fits")
shape, wcs = size(duck), duck.wcs
plot(duck)

##
alm_duck = map2alm(duck)
nm = Pixell.alm2map(alm_duck, size(duck), duck.wcs)
plot(nm)

##
heatmap(duck.data .- nm.data, clim=(-1e-5,1e-5))

##
using PyCall
enmap, curvedsky = pyimport("pixell.enmap"), pyimport("pixell.curvedsky")
imap = pycall(enmap.read_map, PyObject, "test/data/duck.fits")
pyshape = imap.shape
pywcs = imap.wcs
pya = curvedsky.map2alm(imap, lmax=alm_duck.lmax)
newpymap = pycall(enmap.enmap, PyObject, imap,  imap.wcs)
pycall(curvedsky.alm2map, PyObject, pya, newpymap, spin=0)
heatmap(PyArray(imap)' - PyArray(newpymap)', clim=(-1e-5,1e-5))

##
using PowerSpectra

Dl = PowerSpectra.planck_theory_Dl()

##


channelindex(X) = findfirst(first(X), "TEB")

ells = collect(eachindex(Dl[:TT]))
Dlfac = ells .* (ells .+1 ) ./ (2π)

𝐂 = zeros(3, 3, length(ells))
for X in ("T", "E", "B"), Y in ("T", "E", "B")
    c₁, c₂ = channelindex(X), channelindex(Y)
    XY = Symbol(X * Y)
    if XY in keys(Dl)
        𝐂[c₁, c₂ ,:] .= parent(Dl[XY]) ./ Dlfac
        𝐂[c₁, c₂, 1:2] .= 𝐂[c₁, c₂,3]
        𝐂[c₁, c₂, end-10:end] .= 𝐂[c₁, c₂,end-11]
    end
end

lmax = 720
alms = [Alm(lmax, lmax) for s in 1:3]
synalm!(𝐂[:,:,1:lmax+1], alms)

##
m = alm2map(alms[1], shape, wcs)
plot(m)
# plot(abs.(alm2cl(alms[1], alms[2])), yscale=:log10)
# plot!(abs.(𝐂[1,2,1:lmax]))


##


##