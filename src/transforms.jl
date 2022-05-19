
"""Number of pixels in full CAR ring of this WCS."""
fullringsize(wcs::AbstractWCSTransform) = 
    round(Int64, abs(2π / (getunit(wcs) * getcdelt(wcs)[1])))

"""Number of rings in a fullsky version of this WCS."""
fullringnum(wcs::AbstractWCSTransform) = 
    1 + round(Int64, abs(π / (getunit(wcs) * getcdelt(wcs)[2])))

"""Get index of ring given angle"""
function first_last_rings_in_fullsky(shape, wcs::AbstractWCSTransform)
    Δδ = getcdelt(wcs)[end] * getunit(wcs)
    Δθ = abs(Δδ)
    δ₁ = pix2sky(shape, wcs, 1, 1)[2]
    δ₂ = pix2sky(shape, wcs, 1, shape[2])[2]
    θ₁ = π/2 - δ₁
    θ₂ = π/2 - δ₂
    i1 = round(Int64, θ₁ / Δθ) + 1
    i2 = round(Int64, θ₂ / Δθ) + 1

    return i1:sign(i2-i1):i2
end

# libsharp wants ascending colatitude (θ) and increasing right ascension (ϕ)
function get_flip_slices(shape, wcs)
    Δα, Δδ = getcdelt(wcs) .* getunit(wcs)
    flipx_slice = (Δα ≥ 0) ? (1:shape[1]) : (shape[1]:-1:1)
    flipy_slice = (Δδ ≤ 0) ? (1:shape[2]) : (shape[2]:-1:1)
    return flipx_slice, flipy_slice 
end

"""Clenshaw-Curtis quadrature weights, make a Libsharp GeomInfo"""
function make_cc_geom_info(shape, wcs₀::AbstractWCSTransform)

    flipx_slice, flipy_slice = get_flip_slices(shape, wcs₀)
    _, wcs = slice_geometry(shape, wcs₀, flipx_slice, flipy_slice, shape[3:end])
    subinds = first_last_rings_in_fullsky(shape, wcs)
    @assert first(subinds) ≤ last(subinds)  # vertical angle must be increasing
    nringstot = fullringnum(wcs)
    ppring = fullringsize(wcs)
    phi0 = pix2sky(shape, wcs, 1, 2)[1]  # since RA = ϕ

    geom_info_ptr = Ref{Ptr{Cvoid}}()
    μ = chebyshevjacobimoments1(Float64, nringstot, 0.0, 0.0)
    weights = Cdouble.(clenshawcurtisweights(μ)[subinds] .* 2π ./ ppring)
    theta = Cdouble.(range(0, π, length=nringstot)[subinds])

    # constant stuff
    nsubrings = length(subinds)
    phi0s = fill(Cdouble(phi0), nsubrings)
    nph = fill(Cint(ppring), nsubrings)
    stride = ones(Cint, nsubrings)
    offsets = Cptrdiff_t.(ppring .* (0:(nsubrings-1)))

    ccall(
        (:sharp_make_geom_info, Libsharp.libsharp2),
        Cvoid,
        (Cint,       Ref{Cint}, Ref{Cptrdiff_t}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Ptr{Cvoid}}),
        Cint(nsubrings),nph,  offsets,           stride,    phi0s,         theta,        weights,      geom_info_ptr
    )

    return GeomInfo(geom_info_ptr[])
end

"""Complete rings for spherical harmonic transforms. Also converts to AbstractArray{Float64,2}."""
function create_sht_band(m::Enmap)
    shape = size(m)
    flipx_slice, flipy_slice = get_flip_slices(shape, getwcs(m))
    _, wcs = slice_geometry(shape, getwcs(m), flipx_slice, flipy_slice, shape[3:end])
    target_size = (fullringsize(wcs), size(m,2), shape[3:end]...)
    converted = Enmap(zeros(Float64, target_size), wcs)
    other_slices = axes(m)[3:end]

    # extend the ring by adding more columns to the end
    converted.data[1:shape[1], :, other_slices...] .= view(m.data, flipx_slice, flipy_slice, other_slices...)
    return converted
end

function create_sht_band(shape, wcs::AbstractWCSTransform)
    target_size = (fullringsize(wcs), shape[2:end]...)
    return Enmap(zeros(Float64, target_size), wcs)
end

"""Generate an estimate of the Nyquist frequency for a map"""
getlmax(wcs) = fullringsize(wcs) ÷ 2

"""perform forward SHT of an intensity map"""
function map2alm(input_map::Enmap{T,2}; lmax=nothing, mmax=lmax) where T
    if isnothing(lmax)
        lmax = getlmax(getwcs(input_map))
        mmax = lmax
    end
    band = create_sht_band(input_map)
    alm_info = make_triangular_alm_info(lmax, mmax, 1)
    geom_info = make_cc_geom_info(size(band), getwcs(band))
    nalms = alm_count(alm_info)
    npix = map_size(geom_info)
    alm = Alm(lmax, mmax, zeros(ComplexF64, nalms))
    flat_map = reshape(band.data, npix)

    sharp_execute!(
        SHARP_MAP2ALM, 0, 
        [alm.alm], 
        [flat_map],
        geom_info, alm_info, SHARP_DP
    )
    return alm
end


"""perform inverse SHT of an intensity map"""
function alm2map(alm::Alm, shape, wcs::AbstractWCSTransform) where T

    band = create_sht_band(shape, wcs)
    geom_info = make_cc_geom_info(size(band), getwcs(band))
    npix = map_size(geom_info)
    flat_map = reshape(band.data, npix)

    alm_info = make_triangular_alm_info(alm.lmax, alm.mmax, 1)

    sharp_execute!(
        SHARP_ALM2MAP, 0, 
        [alm.alm], 
        [flat_map],
        geom_info, alm_info, SHARP_DP)

    Δα, Δδ = getcdelt(wcs)
    band_width = size(band,1)
    xslice = (Δα ≥ 0) ? (1:shape[1]) : (band_width:-1:(band_width-shape[1]+1))
    yslice = (Δδ ≤ 0) ? (1:shape[2]) : (shape[2]:-1:1)

    return Enmap(band.data[xslice, yslice], wcs)
end

# spin 2 transforms ----

# 2-tuple of scalar enmaps means QU transform
function map2alm(maps::NTuple{2,EM}; lmax=nothing, mmax=lmax) where {T, EM<:Enmap{T,2}}
    if isnothing(lmax)
        lmax = getlmax(getwcs(first(maps)))
        mmax = lmax
    end
    band_maps = create_sht_band.(maps)
    geom_info = make_cc_geom_info(size(band_maps[1]), getwcs(band_maps[1]))
    npix = map_size(geom_info)
    flat_maps = [reshape(m.data, npix) for m in band_maps]

    alm_info = make_triangular_alm_info(lmax, mmax, 1)
    nalm = alm_count(alm_info)
    alms = [Alm(lmax, mmax, zeros(ComplexF64, nalm)) for m in band_maps]

    sharp_execute!(
        SHARP_MAP2ALM, 2, 
        [a.alm for a in alms], 
        flat_maps,
        geom_info, alm_info, SHARP_DP)

    return alms[1], alms[2]
end

# 3-tuple of scalar enmaps means IQU transform
function map2alm(maps::NTuple{3,EM}; lmax=nothing, mmax=lmax) where {T, EM<:Enmap{T,2}}
    if isnothing(lmax)
        lmax = getlmax(getwcs(first(maps)))
        mmax = lmax
    end
    alm_i = map2alm(maps[1]; lmax=lmax, mmax=mmax)
    alm_e, alm_b = map2alm(maps[2:3]; lmax=lmax, mmax=mmax)
    return alm_i, alm_e, alm_b
end

# handle enmaps that are stacks of I, QU, and IQU in the third dim
function map2alm(input_map::Enmap{T,3}; lmax=nothing, mmax=lmax) where T
    if isnothing(lmax)
        lmax = getlmax(getwcs(input_map))
        mmax = lmax
    end
    if size(input_map,3) == 1  # this I map for some reason has a single third dim
        sliced_1d_map = @view input_map[:,:,1]
        return map2alm(sliced_1d_map; lmax=lmax, mmax=mmax)
    elseif size(input_map,3) == 2  # QU map
        return _map2almQU(input_map; lmax=lmax, mmax=mmax)
    elseif size(input_map,3) == 3  # IQU map
        alm_t = map2alm(@view(input_map[:,:,1]); lmax=lmax, mmax=mmax)
        alm_e, alm_b = _map2almQU(@view(input_map[:,:,2:3]); lmax=lmax, mmax=mmax)
        return alm_t, alm_e, alm_b
    end
    throw(ArgumentError("SHTs require shape (nx,ny,ncomp) with 1 ≤ ncomp ≤ 3, for I, QU, and IQU."))
end

# special case QU: have to manually ccall with the columns when 
# Enmap is a stack (nx, ny, 2)
# (internal) perform forward SHT of an QU map
function _map2almQU(input_map::Enmap{T,3}; lmax, mmax) where T

    @assert size(input_map,3) == 2  # check that we have Q and U channels

    band_map = create_sht_band(input_map)
    geom_info = make_cc_geom_info(size(band_map), getwcs(band_map))
    npix = map_size(geom_info)
    flat_map = reshape(band_map.data, npix * 2)

    alm_info = make_triangular_alm_info(lmax, mmax, 1)
    nalm = alm_count(alm_info)
    alm_e = Alm(lmax, mmax, zeros(ComplexF64, nalm))
    alm_b = Alm(lmax, mmax, zeros(ComplexF64, nalm))

    # have to manually ccall here since I didn't add the pointers directly
    GC.@preserve alm_e alm_b flat_map ccall(
        (:sharp_execute, Libsharp.libsharp2),
        Cvoid,
        (Cint, Cint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint,
            Ref{Cdouble}, Ref{Culonglong}),
        SHARP_MAP2ALM, 2, 
        [pointer(alm_e.alm), pointer(alm_b.alm)], 
        [pointer(flat_map, 1), pointer(flat_map, 1+npix)],
        geom_info.ptr, alm_info.ptr, SHARP_DP,
        Ptr{Cdouble}(C_NULL), Ptr{Culonglong}(C_NULL))

    return alm_e, alm_b
end
