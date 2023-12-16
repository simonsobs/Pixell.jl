
abstract type AbstractCARWCS{T} <: AbstractWCSTransform end

"""
Fast custom WCS structure.
"""
struct CarClenshawCurtis{T} <: AbstractCARWCS{T}
    cdelt::Tuple{T,T}
    crpix::Tuple{T,T}
    crval::Tuple{T,T}
    unit::T  # conversion factor to radians
end

struct CarFejer1{T} <: AbstractCARWCS{T}
    cdelt::Tuple{T,T}
    crpix::Tuple{T,T}
    crval::Tuple{T,T}
    unit::T  # conversion factor to radians
end

getunit(T::Type{<:Real}, wcs::AbstractCARWCS) = T(wcs.unit)
getcdelt(wcs::AbstractCARWCS) = wcs.cdelt
getcrpix(wcs::AbstractCARWCS) = wcs.crpix
getcrval(wcs::AbstractCARWCS) = wcs.crval

Base.copy(w::W) where {W<:AbstractCARWCS} = w  # both CAR WCS types are fully immutable

function Base.convert(::Type{CarClenshawCurtis{T}}, w0::WCSTransform) where T
    # FIXME add some asserts
    return CarClenshawCurtis{T}(
        T.(getcdelt(w0)), T.(getcrpix(w0)), T.(getcrval(w0)), getunit(T, w0))
end

function Base.convert(::Type{CarFejer1{T}}, w0::WCSTransform) where T
    # FIXME add some asserts
    return CarFejer1{T}(
        T.(getcdelt(w0)), T.(getcrpix(w0)), T.(getcrval(w0)), getunit(T, w0))
end

function Base.convert(::Type{WCSTransform}, w0::AbstractCARWCS)
    return WCSTransform(2;
        ctype = ["RA---CAR", "DEC--CAR"],
        cdelt = collect(getcdelt(w0)),
        crpix = collect(getcrpix(w0)),
        crval = collect(getcrval(w0)))
end

# this kind of WCS only has two spatial dimensions. this check should be constant-propagated
function Base.getproperty(wcs::AbstractCARWCS, k::Symbol)
    if k == :naxis
        return 2
    end
    return getfield(wcs, k)
end

function Base.show(io::IO, wcs::CarClenshawCurtis{T}) where T
    expr = join(["$k=$(getproperty(wcs, Symbol(k)))"
                 for k in ["naxis","cdelt","crval","crpix"]], ",")
    print(io, "CarClenshawCurtis{$(T)}($expr)")
end
function Base.show(io::IO, wcs::CarFejer1{T}) where T
    expr = join(["$k=$(getproperty(wcs, Symbol(k)))"
                 for k in ["naxis","cdelt","crval","crpix"]], ",")
    print(io, "CarFejer1{$(T)}($expr)")
end


"""
    pix2sky!(m::Enmap, pixcoords, skycoords)

Convert 1-indexed pixels to sky coordinates, in-place. While typically WCS can specify custom units, Pixell 
uses radians throughout.

# Arguments:
- `m::Enmap`: the map that provides a coordinate system
- `pixcoords`: pixel coordinates should be a 2-d array where "pixcoords[:, i]" is the i-th 
    set of coordinates, or a 1-d array representing a single set of coordinates. 
- `skycoords`: output array for sky coordinates, must be same same as pixcoords

# Returns: 
- `Array`: skycoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
       pixcoords =  100 .* rand(2,4096 * 2)
       skycoords =  similar(pixcoords)

julia> pix2sky!(shape, wcs, pixcoords, skycoords)
```
"""
function pix2sky!(shape, wcs::AbstractCARWCS, pixcoords::AbstractArray{TP,2}, 
                  skycoords::AbstractArray{TS,2}; safe=true) where {TP,TS}

    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)

    # compute RA (α) and DEC (δ)
    @inbounds for ipix ∈ axes(pixcoords, 2)
        iα = pixcoords[1, ipix]
        iδ = pixcoords[2, ipix]
        α = α₀ + (iα - iα₀) * Δα
        δ = δ₀ + (iδ - iδ₀) * Δδ
        skycoords[1, ipix] = α
        skycoords[2, ipix] = δ
    end
    
    if safe
        unwind!(skycoords; dims=2)
    end

    return skycoords
end

# the not-in-place version just creates an output array and calls the in-place one above
function pix2sky(shape, wcs::AbstractCARWCS, 
                 pixcoords::AbstractArray{TP,2}; safe=true) where TP
    skycoords = similar(pixcoords)
    return pix2sky!(shape, wcs, pixcoords, skycoords; safe=safe)
end

"""
    pix2sky(m::Enmap{T,N,AA,<:AbstractCARWCS}, ra_pixel, dec_pixel)

Compute the sky position of a single position on the sky.

Only implemented for CAR projections, so
the input map is of type `Enmap{T,N,AA,<:AbstractCARWCS}`.
This takes pixel indices for RA and DEC, and returns a tuple containing
the corresponding RA and DEC.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> pix2sky(shape, wcs, 30.0, 80.0)
(151.0, -11.0)
```
"""
function pix2sky(shape, wcs::AbstractCARWCS, ra_pixel, dec_pixel; safe=true)
    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)
    α = α₀ .+ (ra_pixel .- iα₀) .* Δα
    δ = δ₀ .+ (dec_pixel .- iδ₀) .* Δδ
    if safe
        return rewind(α), rewind(δ)  # FIXME to unwind
    end
    return α, δ
end

# when passing a length-2 vector [ra, dec], return a vector. wraps the pix2sky(shape, wcs, ra_pix, dec_pix)
function pix2sky(shape, wcs::AbstractCARWCS, pixcoords::AbstractVector; safe=true)
    @assert length(pixcoords) == 2
    skycoords = collect(pix2sky(shape, wcs, first(pixcoords), last(pixcoords)))
    if safe
        unwind!(skycoords; dims=2)
    end
    return skycoords
end


function sky2pix!(shape, wcs::AbstractCARWCS, skycoords::AbstractArray{TS,2}, 
                  pixcoords::AbstractArray{TP,2}; safe=true) where {TS,TP}

    # retrieve WCS info
    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)
    Δα⁻¹, Δδ⁻¹ = 1 / Δα, 1 / Δδ

    # compute RA (α) index and DEC (δ) index
    @inbounds for ipix ∈ axes(pixcoords, 2)
        α = skycoords[1, ipix]
        δ = skycoords[2, ipix]
        iα = iα₀ + (α - α₀) * Δα⁻¹
        iδ = iδ₀ + (δ - δ₀) * Δδ⁻¹
        pixcoords[1, ipix] = iα
        pixcoords[2, ipix] = iδ
    end

    if safe
        center_pix = shape[begin:begin+1] ./ 2 .+ 1
        pix_periods = abs.(2π ./ (Δα, Δδ))
        rewind!(view(pixcoords, 1, :); period=pix_periods[1], ref_angle=center_pix[1])
        rewind!(view(pixcoords, 2, :); period=pix_periods[2], ref_angle=center_pix[2])
    end

    return pixcoords
end

# the not-in-place version just creates an output array and calls the in-place one above
function sky2pix(shape, wcs::AbstractCARWCS, 
                 skycoords::AbstractArray{TS,2}; safe=true) where {TS}
    pixcoords = similar(skycoords)
    return sky2pix!(shape, wcs, skycoords, pixcoords; safe=safe)
end


"""
    sky2pix(m::Enmap{T,N,AA,<:AbstractCARWCS}, ra, dec)

Compute 1-indexed pixels into sky coordinates.

Only implemented for CAR (Clenshaw-Curtis variant) projections. Takes 
RA and DEC and returns a tuple containing the corresponding pixel
indices. If vectors of RA and DEC are given, then vectors of 
pixel indices will be returned.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> sky2pix(shape, wcs, deg2rad(30.0), deg2rad(80.0))
(151.0, 171.0)
```
"""
function sky2pix(shape, wcs::AbstractCARWCS, ra::Number, dec::Number; safe=true)
    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)
    pix_ra = iα₀ + (ra - α₀) / Δα
    pix_dec = iδ₀ + (dec - δ₀) / Δδ

    if safe
        center_pix = shape[begin:begin+1] ./ 2 .+ 1
        pix_ra = rewind(pix_ra; period=abs(2π / Δα), ref_angle=center_pix[1])
        pix_dec = rewind(pix_dec; period=abs(2π / Δδ), ref_angle=center_pix[2])
    end
    return pix_ra, pix_dec
end
function sky2pix(shape, wcs::AbstractCARWCS, 
                 ra::AV, dec::AV; safe=true) where {AV<:AbstractVector}
    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)
    Δα⁻¹, Δδ⁻¹ = 1 / Δα, 1 / Δδ
    pix_ra = iα₀ .+ (ra .- α₀) .* Δα⁻¹
    pix_dec = iδ₀ .+ (dec .- δ₀) .* Δδ⁻¹
    
    if safe
        center_pix = shape[begin:begin+1] ./ 2 .+ 1
        rewind!(pix_ra; period=abs(2π * Δα⁻¹), ref_angle=center_pix[1])
        rewind!(pix_dec; period=abs(2π * Δδ⁻¹), ref_angle=center_pix[2])
    end

    return pix_ra, pix_dec
end

# when passing a vector [ra, dec], return a vector. wraps the sky2pix(shape, wcs, ra, dec).
function sky2pix(shape, wcs::AbstractCARWCS, skycoords::AbstractVector; safe=true)
    @assert length(skycoords) == 2
    return collect(sky2pix(shape, wcs::AbstractCARWCS, 
        first(skycoords), last(skycoords); safe=safe))
end

skyarea(shape, wcs::AbstractCARWCS) = skyarea_cyl(shape, wcs)
iscyl(wcs::AbstractCARWCS) = true

"""Generate a similar Enmap whose pixel values are the areas of the pixels in steradians."""
function pixareamap(m::Enmap{T,N,AA,W}) where {T,N,AA,W<:AbstractCARWCS}
    pixareas = similar(m)
    pixareamap!(pixareas)
end

function pixareamap(shape, wcs::AbstractCARWCS) 
    pixareas = Enmap(Array{Float64}(undef, shape), wcs)
    pixareamap!(pixareas)
end

function sliced_wcs(wcs::W, cdelt′, crpix′) where {T, W<:AbstractCARWCS{T}}
    new_wcs = W(cdelt′, crpix′, wcs.crval, wcs.unit)
    return new_wcs
end

function pad(m::Enmap, npix_ra::Int, npix_dec::Int; mode=:center)
    if mode == :center
        center_pad(m, npix_ra, npix_dec)
    elseif mode == :corner
        corner_pad(m, npix_ra, npix_dec)
    end
end

function center_pad(m::Enmap{T,N,AA,W}, npix_ra::Int, npix_dec::Int) where {T,N,AA,W<:AbstractCARWCS}
    new_shape, new_wcs = center_pad(size(m), m.wcs, npix_ra, npix_dec)
    arr = AA(undef, new_shape)
    fill!(arr, zero(T))
    @views arr[
        (begin+npix_ra):(end-npix_ra),
        (begin+npix_dec):(end-npix_dec),
        map(Base.oneto, size(m)[(begin+2):end])...] .= m.data
    return Enmap(arr, new_wcs)
end

function center_pad(shape, wcs::W, npix_ra::Int, npix_dec) where {T,W<:AbstractCARWCS{T}}
    new_shape = (shape[1] + 2npix_ra, shape[2] + 2npix_dec, shape[3:end]...)
    new_wcs = W(
        wcs.cdelt, 
        wcs.crpix .+ (npix_ra, npix_dec), 
        wcs.crval, 
        wcs.unit)  # degree conversion
    return new_shape, new_wcs
end


function corner_pad(m::Enmap{T,N,AA,W}, npix_ra::Int, npix_dec::Int) where {T,N,AA,W<:AbstractCARWCS}
    new_shape, new_wcs = corner_pad(size(m), m.wcs, npix_ra, npix_dec)
    arr = AA(undef, new_shape)
    fill!(arr, zero(T))
    @views arr[
        begin:(end-npix_ra),
        begin:(end-npix_dec),
        map(Base.oneto, size(m)[(begin+2):end])...] .= m.data
    return Enmap(arr, new_wcs)
end

function corner_pad(shape, wcs::AbstractCARWCS{T}, npix_ra::Int, npix_dec) where T
    new_shape = (shape[1] + npix_ra, shape[2] + npix_dec, shape[3:end]...)
    new_wcs = deepcopy(wcs)
    return new_shape, new_wcs
end

pad(m::Enmap{T,N,AA,W}, npix; kwargs...) where {T,N,AA,W<:AbstractCARWCS} = pad(m, npix, npix; kwargs...)

