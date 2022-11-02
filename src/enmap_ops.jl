
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

Convert 1-indexed pixels to sky coordinates. While typically WCS can specify custom units, Pixell 
uses radians throughout.

# Arguments:
- `m::Enmap`: the map that provides a coordinate system
- `pixcoords`: pixcoords should be a 2-d array where "pixcoords[:, i]" is the i-th set of coordinates, 
    or a 1-d array representing a single set of coordinates. 

# Returns: 
- `Array`: same shape as pixcoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
julia> pix2sky(shape, wcs, [1.0, 1.0])
2-element Vector{Float64}:
 -3.141592653589793
 -1.5707963267948966
```
"""
function pix2sky end

pix2sky(m::Enmap, p1; safe=true) = pix2sky(size(m), getwcs(m), p1; safe=safe)
pix2sky(m::Enmap, p1, p2; safe=true) = pix2sky(size(m), getwcs(m), p1, p2; safe=safe)
sky2pix(m::Enmap, p1; safe=true) = sky2pix(size(m), getwcs(m), p1; safe=safe)
sky2pix(m::Enmap, p1, p2; safe=true) = sky2pix(size(m), getwcs(m), p1, p2; safe=safe)

pix2sky!(m::Enmap, p1, p2; safe=true) = pix2sky!(size(m), getwcs(m), p1, p2; safe=safe)
sky2pix!(m::Enmap, p1, p2; safe=true) = sky2pix!(size(m), getwcs(m), p1, p2; safe=safe)



"""
    sky2pix(m::Enmap, skycoords)

Convert sky coordinates to 1-indexed pixels. While typically WCS can specify custom units, Pixell 
uses radians throughout.

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

Convert sky coordinates to 1-indexed pixels, in-place. While typically WCS can specify custom units, Pixell 
uses radians throughout.

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
function sky2pix! end




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


function pad(m::Enmap{T,N,AA}, npix::Int) where {T,N,AA}
    new_shape, new_wcs = pad(size(m), m.wcs, npix)
    arr = AA(undef, new_shape)
    return Enmap(arr, new_wcs)
end

struct SkyBoundingBox{T}
    α_min::T
    δ_min::T
    α_max::T
    δ_max::T
end

function SkyBoundingBox(c1::Tuple{T,T}, c2::Tuple{T,T}) where T
    α_min, α_max = min(c1[1], c2[1]), max(c1[1], c2[1])
    δ_min, δ_max = min(c1[2], c2[2]), max(c1[2], c2[2])
    return SkyBoundingBox{T}(α_min, δ_min, α_max, δ_max)
end

function in(skycoords::Tuple{T,T}, bbox::SkyBoundingBox) where T
    α, δ = skycoords
    b = bbox
    return (b.α_min ≤ α ≤ b.α_max) && (b.δ_min ≤ δ ≤ b.δ_max)
end
