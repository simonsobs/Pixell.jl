
"""
    rewind(angles, period=2π, ref_angle=0)

Given angles or other cyclic coordinates with the specified period, such that
the angle + period has the same meaning as the original angle, this function adds or subtracts 
multiples of the period such that the result has the same meaning, but now all angles lie in
an interval of length the specified period, centered on the reference angle `ref_angle`.
"""
function rewind(angles; period=2π, ref_angle=0)
    half_period = period / 2
    return ref_angle .+ mod.(angles .- ref_angle .+ half_period, period) .- half_period
end

function rewind!(angles; period=2π, ref_angle=0)
    half_period = period / 2
    angles .= ref_angle .+ mod.(angles .- ref_angle .+ half_period, period) .- half_period
    return angles
end

function unwind(angles; dims=nothing, period=2π, ref_angle=0)
    wound_angles = rewind(angles; period=period, ref_angle=ref_angle)
    return unwrap(wound_angles .- ref_angle; dims=dims, range=period) .+ ref_angle
end

function unwind!(angles; dims=nothing, period=2π, ref_angle=0)
    rewind!(angles; period=period, ref_angle=ref_angle)
    angles .-= ref_angle
    unwrap!(angles; dims=dims, range=period)
    angles .+= ref_angle
    return angles
end


"""
    pix2sky(m::Enmap, pixcoords)

Convert 1-indexed pixels to sky coordinates. The output sky coordinates are determined by WCS,
but usually are in units of degrees. 

# Arguments:
- `m::Enmap`: the map that provides a coordinate system
- `pixcoords`: pixcoords should be a 2-d array where "pixcoords[:, i]" is the i-th set of coordinates, 
    or a 1-d array representing a single set of coordinates. 

# Returns: 
- `Array`: same shape as pixcoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> pix2sky(shape, wcs [1.0, 1.0])
2-element Vector{Float64}:
 180.0
 -90.0
```
"""
function pix2sky end

pix2sky(m::Enmap, p1; safe=true) = pix2sky(size(m), getwcs(m), p1; safe=safe)
pix2sky(m::Enmap, p1, p2; safe=true) = pix2sky(size(m), getwcs(m), p1, p2; safe=safe)
sky2pix(m::Enmap, p1; safe=true) = sky2pix(size(m), getwcs(m), p1; safe=safe)
sky2pix(m::Enmap, p1, p2; safe=true) = sky2pix(size(m), getwcs(m), p1, p2; safe=safe)

pix2sky!(m::Enmap, p1, p2; safe=true) = pix2sky!(size(m), getwcs(m), p1, p2; safe=safe)
sky2pix!(m::Enmap, p1, p2; safe=true) = sky2pix!(size(m), getwcs(m), p1, p2; safe=safe)

## The default fallback (no projection specified) is to call WCS, which calls the C library.
## If you wanted to replicate Pixell behavior, add 1 to x and y of pixcoords.
## We implement custom routines for CAR (Clenshaw-Curtis variant) instead of using these.
function pix2sky(shape, wcs::WCSTransform, pixcoords; safe=true)
    angle_unit = getunit(wcs)
    skycoords = pix_to_world(wcs, pixcoords) .* angle_unit
    if safe
        return unwind(skycoords; dims=2)
    end
    return skycoords
end
function pix2sky!(shape, wcs::WCSTransform, pixcoords, skycoords; safe=true)
    angle_unit = getunit(wcs)
    pix_to_world!(wcs, pixcoords, skycoords)
    skycoords .*= angle_unit
    if safe
        unwind!(skycoords; dims=2)
    end
    return skycoords
end
function sky2pix(shape, wcs::WCSTransform, skycoords; safe=true)
    angle_unit = getunit(wcs)
    inverse_angle_unit = 1 / angle_unit
    pixcoords = world_to_pix(wcs, skycoords .* inverse_angle_unit)
    if safe
        center_pix = shape[begin:begin+1] ./ 2 .+ 1
        pix_periods = abs.(2π ./ (getcdelt(wcs) .* angle_unit))
        rewind!(view(pixcoords, 1, :); period=pix_periods[1], ref_angle=center_pix[1])
        rewind!(view(pixcoords, 2, :); period=pix_periods[2], ref_angle=center_pix[2])
    end
    return pixcoords
end
function sky2pix!(shape, wcs::WCSTransform, skycoords, pixcoords; safe=true)
    angle_unit = getunit(wcs)
    inverse_angle_unit = 1 / angle_unit
    world_to_pix!(wcs, skycoords .* inverse_angle_unit, pixcoords)
    if safe
        center_pix = shape[begin:begin+1] ./ 2 .+ 1
        pix_periods = abs.(2π ./ (getcdelt(wcs) .* angle_unit))
        rewind!(view(pixcoords, 1, :); period=pix_periods[1], ref_angle=center_pix[1])
        rewind!(view(pixcoords, 2, :); period=pix_periods[2], ref_angle=center_pix[2])
    end
    return pixcoords
end
pix2sky(shape, wcs::WCSTransform,ra::Number, dec::Number; safe=true) = 
    pix2sky(shape, wcs, Float64[ra, dec]; safe=safe)
sky2pix(shape, wcs::WCSTransform,ra::Number, dec::Number; safe=true) = 
    sky2pix(shape, wcs, Float64[ra, dec]; safe=safe)

function pix2sky(shape, wcs::WCSTransform, ra::AbstractVector, dec::AbstractVector; safe=true)
    result = pix2sky(shape, wcs, convert(Array{Float64,2}, [ra'; dec']); safe=safe)
    return result[1,:], result[2,:]
end
function sky2pix(shape, wcs::WCSTransform, ra::AbstractVector, dec::AbstractVector; safe=true)
    result = sky2pix(shape, wcs, convert(Array{Float64,2}, [ra'; dec']); safe=safe)
    return result[1,:], result[2,:]
end
"""
    pix2sky!(m::Enmap, pixcoords, skycoords)

Convert 1-indexed pixels to sky coordinates, in-place. The output sky coordinates are 
determined by WCS, but usually are in units of degrees. 

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
function pix2sky!(shape, wcs::CarClenshawCurtis, pixcoords::AbstractArray{TP,2}, 
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
function pix2sky(shape, wcs::CarClenshawCurtis, 
                 pixcoords::AbstractArray{TP,2}; safe=true) where TP
    skycoords = similar(pixcoords)
    return pix2sky!(shape, wcs, pixcoords, skycoords; safe=safe)
end

"""
    pix2sky(m::Enmap{T,N,AA,<:CarClenshawCurtis}, ra_pixel, dec_pixel)

Compute the sky position of a single position on the sky.

Only implemented for CAR (Clenshaw-Curtis variant) projections, so
the input map is of type `Enmap{T,N,AA,<:CarClenshawCurtis}`.
This takes pixel indices for RA and DEC, and returns a tuple containing
the corresponding RA and DEC.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> pix2sky(shape, wcs, 30.0, 80.0)
(151.0, -11.0)
```
"""
function pix2sky(shape, wcs::CarClenshawCurtis, ra_pixel, dec_pixel; safe=true)
    angle_unit = getunit(wcs)
    α₀, δ₀ = getcrval(wcs) .* angle_unit
    Δα, Δδ = getcdelt(wcs) .* angle_unit
    iα₀, iδ₀ = getcrpix(wcs)
    α = α₀ .+ (ra_pixel .- iα₀) .* Δα
    δ = δ₀ .+ (dec_pixel .- iδ₀) .* Δδ
    if safe
        return rewind(α), rewind(δ)
    end
    return α, δ
end

# when passing a length-2 vector [ra, dec], return a vector. wraps the pix2sky(shape, wcs, ra_pix, dec_pix)
function pix2sky(shape, wcs::CarClenshawCurtis, pixcoords::AbstractVector; safe=true)
    @assert length(pixcoords) == 2
    skycoords = collect(pix2sky(shape, wcs, first(pixcoords), last(pixcoords)))
    if safe
        unwind!(skycoords; dims=2)
    end
    return skycoords
end


"""
    sky2pix(m::Enmap, skycoords)

Convert sky coordinates to 1-indexed pixels. The input sky coordinates are determined by WCS,
but usually are in units of degrees. 

# Arguments:
- `m::Enmap`: the map to obtain the coordinates from
- `skycoords`: skycoords should be a 2-d array where "skycoords[:, i]" is the i-th set of coordinates, 
    or a 1-d array representing a single set of coordinates. 

# Returns: 
- `Array`: same shape as skycoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> sky2pix(shape, wcs, [30.0, 50.0])
2-element Vector{Float64}:
 151.0
 141.0
```
"""
function sky2pix end

"""
    sky2pix!(m::Enmap, skycoords, pixcoords)

Convert sky coordinates to 1-indexed pixels, in-place. The input sky coordinates are 
determined by WCS, but usually are in units of degrees. 

# Arguments:
- `m::Enmap`: the map that provides a coordinate system
- `skycoords`: sky coordinates should be a 2-d array where "skycoords[:, i]" is the i-th 
    set of coordinates, or a 1-d array representing a single set of coordinates. 
- `pixcoords`: output array for pixel coordinates, must be same same as pixcoords

# Returns: 
- `Array`: pixcoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
       skycoords =  similar(pixcoords)
       pixcoords =  100 .* rand(2,4096 * 2)
julia> sky2pix!(shape, wcs, skycoords, pixcoords)
```
"""
function sky2pix!(shape, wcs::CarClenshawCurtis, skycoords::AbstractArray{TS,2}, 
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
function sky2pix(shape, wcs::CarClenshawCurtis, 
                 skycoords::AbstractArray{TS,2}; safe=true) where {TS}
    pixcoords = similar(skycoords)
    return sky2pix!(shape, wcs, skycoords, pixcoords; safe=safe)
end


"""
    sky2pix(m::Enmap{T,N,AA,<:CarClenshawCurtis}, ra, dec)

Compute 1-indexed pixels into sky coordinates.

Only implemented for CAR (Clenshaw-Curtis variant) projections. Takes 
RA and DEC and returns a tuple containing the corresponding pixel
indices. If vectors of RA and DEC are given, then vectors of 
pixel indices will be returned.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> sky2pix(shape, wcs, 30.0, 80.0)
(151.0, 171.0)
```
"""
function sky2pix(shape, wcs::CarClenshawCurtis, ra::Number, dec::Number; safe=true)
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
function sky2pix(shape, wcs::CarClenshawCurtis, 
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
function sky2pix(shape, wcs::CarClenshawCurtis, skycoords::AbstractVector; safe=true)
    @assert length(skycoords) == 2
    return collect(sky2pix(shape, wcs::CarClenshawCurtis, 
        first(skycoords), last(skycoords); safe=safe))
end

"""Check if a WCS is a cylindrical pixelization"""
function iscyl(wcs::WCSTransform)
    # right now we only support CAR
    ctype = wcs.ctype
    if ctype[1] == "RA---CAR" && ctype[2] == "DEC--CAR"
        return true
    end
    return false
end
iscyl(wcs::CarClenshawCurtis) = true

"""Area of a patch given shape and WCS, in steradians."""
function skyarea(shape, wcs::WCSTransform)
    if iscyl(wcs)
        return skyarea_cyl(shape, wcs)
    end 
    error("Area for this WCS is not implemented.")
end
skyarea(shape, wcs::CarClenshawCurtis) = skyarea_cyl(shape, wcs)

# generic cylindrical area calculation
function skyarea_cyl(shape, wcs::AbstractWCSTransform)
    N_α, N_δ = shape
    α, δ = pix2sky(shape, wcs, SVector(0., 0.), SVector(0.5, N_δ+0.5); safe=false)
    δ₁, δ₂ = sort(δ)
    δ₁, δ₂ = max(-π/2, δ₁), min(π/2, δ₂)
    Δα, Δδ = getcdelt(wcs) .* getunit(wcs)
    return (sin(δ₂) - sin(δ₁)) * abs(Δα) * N_α
end

"""Generate a similar Enmap whose pixel values are the areas of the pixels in steradians."""
function pixareamap(m::Enmap{T,N,AA,W}) where {T,N,AA,W<:CarClenshawCurtis}
    pixareas = similar(m)
    pixareamap!(pixareas)
end

function pixareamap(shape, wcs::CarClenshawCurtis) 
    pixareas = Enmap(Array{Float64}(undef, shape), wcs)
    pixareamap!(pixareas)
end

"""In-place write to pixels the areas of those pixels in steradians."""
function pixareamap!(pixareas::Enmap)
    shape = size(pixareas)
    wcs = getwcs(pixareas)
    ncols, nrows = shape
    Δα, Δδ = abs.(getcdelt(wcs) .* getunit(wcs))

    for i in 1:nrows
        α, δ = pix2sky(shape, wcs, SVector(1., 1.), SVector(i-0.5, i+0.5); safe=false)
        δ₁, δ₂ = sort(δ)
        δ₁, δ₂ = max(-π/2, δ₁), min(π/2, δ₂)
        pixareas[:,i] .= (sin(δ₂) - sin(δ₁)) * Δα
    end

    return pixareas
end




# this set of slices is necessary because colons are somehow not expanded
slice_geometry(shape::Tuple, wcs, sel_x::Colon, s...) = slice_geometry(shape, wcs, 1:shape[1], s...)
slice_geometry(shape::Tuple, wcs, sel_x::Colon, sel_y::Colon, s...) = slice_geometry(shape, wcs, 1:shape[1], 1:shape[2], s...)
slice_geometry(shape::Tuple, wcs, sel_x, sel_y::Colon, s...) = slice_geometry(shape, wcs, sel_x, 1:shape[2], s...)

# this set of slices is neccesary because slices that remove one of the first two dimensions are not allowed, because
# our WCS structure ALWAYS has two axes. However, it can be desirable to write code like "enmap[1,:] .= 0.0" instead 
# even if the object on the left is no longer an Enmap
slice_geometry(shape::Tuple, wcs, sel_x::Integer, sel_y::UnitRange, s...) = slice_geometry(shape, wcs, sel_x:sel_x, sel_y, s...)
slice_geometry(shape::Tuple, wcs, sel_x::Integer, sel_y::Integer, s...) = slice_geometry(shape, wcs, sel_x:sel_x, sel_y:sel_y, s...)
slice_geometry(shape::Tuple, wcs, sel_x::UnitRange, sel_y::Integer, s...) = slice_geometry(shape, wcs, sel_x, sel_y:sel_y, s...)

# does not implement "nowrap" like pixell
function slice_geometry(shape_all::Tuple, wcs, sel_x, sel_y, sel_others...)
    other_dims = shape_all[3:end]
    sel = (sel_x, sel_y)

    starts = map(s -> (step(s) > 0) ? first(s) - 1 : first(s), sel)  # handle backward ranges too
    steps = step.(sel)
    sel_sizes = last.(sel) .- first.(sel) .+ steps
    crpix′ = (getcrpix(wcs) .- (starts .+ 0.5)) ./ steps .+ 0.5
    cdelt′ = getcdelt(wcs) .* steps
    shape = sel_sizes .÷ steps
    wcs′ = sliced_wcs(wcs, cdelt′, crpix′)

    return (shape..., other_dims...), wcs′
end
slice_geometry(m::Enmap, sel_all::Vararg) = slice_geometry(size(m), getwcs(m), sel_all...)


function sliced_wcs(wcs::WCSTransform, cdelt′, crpix′)
    new_wcs = deepcopy(wcs)
    new_wcs.cdelt = collect(cdelt′)
    new_wcs.crpix = collect(crpix′)
    return new_wcs
end

function sliced_wcs(wcs::CarClenshawCurtis{T}, cdelt′, crpix′) where T
    new_wcs = CarClenshawCurtis{T}(cdelt′, crpix′, wcs.crval, wcs.unit)
    return new_wcs
end
