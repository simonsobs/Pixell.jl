
"""Number of pixels in full CAR ring of this WCS."""
fullringsize(wcs::AbstractWCSTransform) = 
    round(Int64, abs(2π / (get_unit(wcs) * cdelt(wcs)[1])))

"""Number of rings in a fullsky version of this WCS."""
fullringnum(wcs::AbstractWCSTransform) = 
    1 + round(Int64, abs(π / (get_unit(wcs) * cdelt(wcs)[2])))

"""Get index of ring given angle"""
function first_last_rings_in_fullsky(shape, wcs::AbstractWCSTransform)
    Δδ = cdelt(wcs)[end] * get_unit(wcs)
    Δθ = abs(Δδ)
    δ₁ = pix2sky(shape, wcs, 1, 1)[end]
    δ₂ = pix2sky(shape, wcs, 1, shape[2])[end]
    θ₁ = π/2 - δ₁
    θ₂ = π/2 - δ₂
    i1 = round(Int64, θ₁ / Δθ) + 1
    i2 = round(Int64, θ₂ / Δθ) + 1
    return i1:i2
end

function get_flip_slices(shape, wcs)
    Δα, Δδ = cdelt(wcs) .* get_unit(wcs)
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

    GeomInfo(geom_info_ptr[])
end


"""Complete rings for spherical harmonic transforms. Also converts to AbstractArray{Float64,2}."""
function create_sht_band(m::Enmap)
    shape = size(m)
    flipx_slice, flipy_slice = get_flip_slices(shape, getwcs(m))
    _, wcs = slice_geometry(shape, getwcs(m), flipx_slice, flipy_slice, shape[3:end])
    target_size = (fullringsize(wcs), size(m,2))
    converted = Enmap(zeros(Float64, target_size), wcs)
    converted.data[1:shape[1],:] .= view(m.data, flipx_slice, flipy_slice)
    return converted
end


getlmax(wcs) = 3 * fullringsize(wcs)

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
    map1d = reshape(band.data, npix)

    sharp_execute!(
        SHARP_MAP2ALM, 0, 
        [alm.alm], 
        [map1d],
        geom_info, alm_info, SHARP_DP
    )
    return alm
end

