
function create_car_wcs(::Type{WCSTransform}, cdelt, crpix, crval)
    wcs = WCSTransform(2;
        ctype = ["RA---CAR", "DEC--CAR"],
        cdelt = collect(cdelt),
        crpix = collect(crpix),
        crval = collect(crval))

    return wcs
end

# try to follow WCS.WCSTransform with degree units
function create_car_wcs(::Type{CAR{T}}, cdelt, crpix, crval) where T
    return CAR{T}(
        (cdelt[1], cdelt[2]), 
        (crpix[1], crpix[2]), 
        (crval[1], crval[2]), 
        π/180)  # degree conversion
end
create_car_wcs(::Type{CAR}, cdelt, crpix, crval) = 
    create_car_wcs(CAR{Float64}, cdelt, crpix, crval)

"""
    fullsky_geometry([W=CAR], res; shape = nothing, dims = ())

Generates a full-sky geometry.

# Arguments:
- `proj=CAR()`: [optional] projection
- `res`: resolution in radians. Passing a Number produces a square pixel.
Passing a tuple with (ΔRA, ΔDEC) produces a rectangular pixel.

# Keywords
- `shape::NTuple=nothing`: shape of the map. If not specified, will be computed.
- `dims::NTuple=()`: additional dimensions to append to the shape, such as (3,) for IQU
to generate a map with `(nx, ny, 3)`.

# Returns: 
- `shape::Tuple, wcs::W`: a tuple containing the shape of the map and the WCS

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(30/60))  # 30 arcmin pixel
((720, 361), WCSTransform(naxis=2,cdelt=[-0.5, 0.5],crval=[0.25, 0.0],crpix=[360.5, 181.0]))
```
"""
function fullsky_geometry(W::Type{<:AbstractWCSTransform}, res; shape = nothing, dims = ())
    if isnothing(shape)
        shape = (round.(Int, (2π, π) ./ res .+ (0, 1)))  # CAR has pixels on poles
    end
    nx, ny = shape
    res_in_radians = ustrip.(uconvert.(radian, res))
    resx, resy = res_in_radians

    @assert abs(resx * nx - 2π) < 1e-8 "Horizontal resolution does not evenly divide the sky; this is required for SHTs."
    @assert abs(resy * (ny - 1) - π) < 1e-8 "Vertical resolution does not evenly divide the sky; this is required for SHTs."

    # Note the reference point is shifted by half a pixel to keep
    # the grid in bounds, from ra=180+cdelt/2 to ra=-180+cdelt/2.
    wcs = create_car_wcs(W,  
        (-360.0 / nx, 180.0 / (ny - 1)),      # cdelt
        (floor(nx / 2) + 0.5, (ny + 1) / 2),  # crpix
        (resy * 90 / π, 0.0)                  # crval
    )

    return (nx, ny, dims...), wcs
end

fullsky_geometry(W::Type{<:AbstractWCSTransform}, res::Number; shape = nothing, dims = ()) =
    fullsky_geometry(W, (res, res); shape = shape, dims = dims)

fullsky_geometry(res; shape = nothing, dims = ()) =
    fullsky_geometry(CAR{Float64}, res; shape = shape, dims = dims)


# ONLY DOES CAR FOR NOW
function geometry(W::Type{<:AbstractWCSTransform}, bbox_coords, res)

    # get everything into radians
    res_in_radians = ustrip.(uconvert.(radian, res))
    resx, resy = res_in_radians
    @assert abs(2π / resx - round(2π / resx)) < 1e-8 "Horizontal resolution" * 
        " does not evenly divide the sky; this is required for SHTs."
    @assert abs(2π / resy - round(2π / resy)) < 1e-8 "Vertical resolution" *
        " does not evenly divide the sky; this is required for SHTs."

    # unpack the bounding box coordinates and set up the shape
    pos1 = ustrip.(uconvert.(radian, bbox_coords[:,1]))
    pos2 = ustrip.(uconvert.(radian, bbox_coords[:,2]))
    Δra, Δdec = abs.(pos1 .- pos2)
    shape = round.(Int, (Δra, Δdec) ./ res_in_radians)

    # set up coordinates such that the [1, 1] coordinate is the first pos of bbox
    mid = (pos1 .+ pos2) ./ 2
    crval = (mid[1], 0.0)
    cdelt = abs.(res_in_radians) .* sign.(pos2 .- pos1)
    crpix = 1 .- (pos1 .- crval) ./ cdelt

    wcs = create_car_wcs(W,  
        rad2deg.(cdelt),  # cdelt
        crpix,            # crpix
        rad2deg.(crval))  # crval

    return shape, wcs
end

geometry(W::Type{<:AbstractWCSTransform}, bbox_coords, res::Number) = 
    geometry(W, bbox_coords, (res, res))
