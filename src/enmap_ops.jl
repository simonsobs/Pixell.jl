

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
       m = Enmap(rand(shape...), wcs)
julia> pix2sky(m, [1.0, 1.0])
2-element Vector{Float64}:
 180.0
 -90.0
```
"""
function pix2sky end

## The default fallback (no projection specified) is to call WCS, which calls the C library.
## If you wanted to replicate Pixell behavior, add 1 to x and y of pixcoords.
## We implement custom routines for CAR (Clenshaw-Curtis variant) instead of using these.
function pix2sky(m::Enmap{T}, pixcoords) where T
    angle_unit = get_unit(T, w)
    return pix_to_world(getwcs(m), pixcoords) .* angle_unit
end
function pix2sky!(m::Enmap{T}, pixcoords, skycoords) where T
    angle_unit = get_unit(T, w)
    pix_to_world!(getwcs(m), pixcoords, skycoords)
    skycoords .*= angle_unit
    return skycoords
end
function sky2pix(m::Enmap{T}, skycoords) where T
    inverse_angle_unit = 1 / get_unit(T, w)
    return world_to_pix(getwcs(m), skycoords .* inverse_angle_unit)
end
function sky2pix!(m::Enmap{T}, skycoords, pixcoords) where T
    inverse_angle_unit = 1 / get_unit(T, w)
    return world_to_pix!(getwcs(m), skycoords .* inverse_angle_unit, pixcoords)
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
       m = Enmap(rand(shape...), wcs)
       pixcoords =  100 .* rand(2,4096 * 2)
       skycoords =  similar(pixcoords)

julia> pix2sky!(m, pixcoords, skycoords)
```
"""
function pix2sky!(m::Enmap{T,N,AA,CarClenshawCurtis},
    pixcoords::AbstractArray{TP,2}, skycoords::AbstractArray{TS,2}
) where {T,N,AA<:AbstractArray{T,N},TP,TS}

    # retrieve WCS info
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)

    # compute RA (α) and DEC (δ)
    @inbounds for ipix ∈ axes(pixcoords, 2)
        iα = pixcoords[1, ipix]
        iδ = pixcoords[2, ipix]
        α = α₀ + (iα - iα₀) * Δα
        δ = δ₀ + (iδ - iδ₀) * Δδ
        skycoords[1, ipix] = α
        skycoords[2, ipix] = δ
    end

    return skycoords
end

# the not-in-place version just creates an output array and calls the in-place one above
function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, pixcoords::AbstractArray{TP,2}
) where {T,N,AA<:AbstractArray{T,N},TP}
    skycoords = similar(pixcoords)
    return pix2sky!(m, pixcoords, skycoords)
end

"""
    pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, ra_pixel, dec_pixel)

Compute the sky position of a single position on the sky.

Only implemented for CAR (Clenshaw-Curtis variant) projections, so
the input map is of type `Enmap{T,N,AA,CarClenshawCurtis}`.
This takes pixel indices for RA and DEC, and returns a tuple containing
the corresponding RA and DEC.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
       m = Enmap(rand(shape...), wcs)
julia> pix2sky(m, 30.0, 80.0)
(151.0, -11.0)
```
"""
function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 ra_pixel::Number, dec_pixel::Number) where {T,N,AA<:AbstractArray{T,N}}
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)
    α = α₀ + (ra_pixel - iα₀) * Δα
    δ = δ₀ + (dec_pixel - iδ₀) * Δδ
    return α, δ
end
function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 ra_pixel::AV, dec_pixel::AV) where {T,N,AA<:AbstractArray{T,N}, AV<:AbstractVector}
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)
    α = α₀ .+ (ra_pixel .- iα₀) .* Δα
    δ = δ₀ .+ (dec_pixel .- iδ₀) .* Δδ
    return α, δ
end

# when passing a length-2 vector [ra, dec], return a vector. wraps the pix2sky(m, ra_pix, dec_pix)
function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 pixcoords::AbstractVector) where {T,N,AA<:AbstractArray{T,N}}
    @assert length(pixcoords) == 2
    return collect(pix2sky(m, first(pixcoords), last(pixcoords)))
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
       m = Enmap(rand(shape...), wcs)
julia> sky2pix(m, [30.0, 50.0])
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
       m = Enmap(rand(shape...), wcs)
       skycoords =  similar(pixcoords)
       pixcoords =  100 .* rand(2,4096 * 2)
julia> sky2pix!(m, skycoords, pixcoords)
```
"""
function sky2pix!(m::Enmap{T,N,AA,CarClenshawCurtis}, skycoords::AbstractArray{TS,2}, 
                  pixcoords::AbstractArray{TP,2}) where {T,N,AA<:AbstractArray{T,N},TS,TP}

    # retrieve WCS info
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)
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

    return pixcoords
end

# the not-in-place version just creates an output array and calls the in-place one above
function sky2pix(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 skycoords::AbstractArray{TS,2}) where {T,N,AA<:AbstractArray{T,N},TS}
    pixcoords = similar(skycoords)
    return sky2pix!(m, skycoords, pixcoords)
end


"""
    sky2pix(m::Enmap{T,N,AA,CarClenshawCurtis}, ra, dec)

Compute 1-indexed pixels into sky coordinates.

Only implemented for CAR (Clenshaw-Curtis variant) projections. Takes 
RA and DEC and returns a tuple containing the corresponding pixel
indices. If vectors of RA and DEC are given, then vectors of 
pixel indices will be returned.

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
       m = Enmap(rand(shape...), wcs)
julia> sky2pix(m, 30.0, 80.0)
(151.0, 171.0)
```
"""
function sky2pix(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 ra::Number, dec::Number) where {T,N,AA<:AbstractArray{T,N}}
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)
    pix_ra = iα₀ + (ra - α₀) / Δα
    pix_dec = iδ₀ + (dec - δ₀) / Δδ
    return pix_ra, pix_dec
end
function sky2pix(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 ra::AV, dec::AV) where {T,N,AA<:AbstractArray{T,N}, AV<:AbstractVector}
    w = getwcs(m)
    angle_unit = get_unit(T, w)
    α₀, δ₀ = crval(w) .* angle_unit
    Δα, Δδ = cdelt(w) .* angle_unit
    iα₀, iδ₀ = crpix(w)
    Δα⁻¹, Δδ⁻¹ = 1 / Δα, 1 / Δδ
    pix_ra = iα₀ .+ (ra .- α₀) .* Δα⁻¹
    pix_dec = iδ₀ .+ (dec .- δ₀) .* Δδ⁻¹
    return pix_ra, pix_dec
end

# when passing a vector [ra, dec], return a vector. wraps the sky2pix(m, ra, dec).
function sky2pix(m::Enmap{T,N,AA,CarClenshawCurtis},
                 skycoords::AbstractVector) where {T,N,AA<:AbstractArray{T,N}}
    @assert length(skycoords) == 2
    return collect(sky2pix(m, first(skycoords), last(skycoords)))
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
    crpix′ = (crpix(wcs) .- (starts .+ 0.5)) ./ steps .+ 0.5
    cdelt′ = cdelt(wcs) .* steps
    shape = sel_sizes .÷ steps

    wcs′ = deepcopy(wcs)
    wcs′.cdelt = collect(cdelt′)
    wcs′.crpix = collect(crpix′)
    return (shape..., other_dims...), wcs′
end
slice_geometry(m::Enmap, sel_all::Vararg) = slice_geometry(size(m), getwcs(m), sel_all...)
