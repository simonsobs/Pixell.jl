
"""
fullsky_geometry([P=CarClenshawCurtis], res; shape = nothing, dims = ())

Generates a full-sky geometry.

# Arguments:
- `proj=CarClenshawCurtis`: [optional] projection
- `res`: resolution in radians. Passing a Number produces a square pixel.
Passing a tuple with (ΔRA, ΔDEC) produces a rectangular pixel.

# Keywords
- `shape::NTuple=nothing`: shape of the map. If not specified, will be computed.
- `dims::NTuple=()`: additional dimensions to append to the shape, such as (3,) for IQU
to generate a map with `(nx, ny, 3)`.

# Returns: 
- `shape::Tuple, wcs::WCSTransform`: a tuple containing the shape of the map and the WCS

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(30/60))  # 30 arcmin pixel
((720, 361), WCSTransform(naxis=2,cdelt=[-0.5, 0.5],crval=[0.25, 0.0],crpix=[360.5, 181.0]))
```
"""
function fullsky_geometry(::Type{<:CarClenshawCurtis}, res; shape = nothing, dims = ())
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
    wcs = WCSTransform(2;
        cdelt = [-360.0 / nx, 180.0 / (ny - 1)],
        ctype = ["RA---CAR", "DEC--CAR"],
        crpix = [floor(nx / 2) + 0.5, (ny + 1) / 2],
        crval = [resy * 90 / π, 0])

    return (nx, ny, dims...), wcs
end

fullsky_geometry(proj::Type{<:CarClenshawCurtis}, res::Number; shape = nothing, dims = ()) =
    fullsky_geometry(proj, (res, res); shape = shape, dims = dims)

fullsky_geometry(res; shape = nothing, dims = ()) =
    fullsky_geometry(CarClenshawCurtis, res; shape = shape, dims = dims)


function geometry(::Type{<:CarClenshawCurtis}, bbox_coords, res)

    # get everything into radians
    res_in_radians = ustrip.(uconvert.(radian, res))
    resx, resy = res_in_radians
    @show res_in_radians resx resy
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
    crval = [mid[1], 0.0]
    cdelt = abs.(res_in_radians) .* sign.(pos2 .- pos1)
    crpix = 1 .- (pos1 .- crval) ./ cdelt

    # WCS is in degrees by default, now we have to convert every angle
    wcs = WCSTransform(2;
        cdelt = rad2deg.(cdelt),
        ctype = ["RA---CAR", "DEC--CAR"],
        crpix = crpix,
        crval = rad2deg.(crval))

    return shape, wcs
end

geometry(proj::Type{<:CarClenshawCurtis}, bbox_coords, res::Number) = 
    geometry(proj, bbox_coords, (res, res))
