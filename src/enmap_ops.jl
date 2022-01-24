
function fullsky_geometry(proj::Type{<:CarClenshawCurtis}, res; shape = nothing, dims = ())
    if isnothing(shape)
        shape = (round.(Int, (2π, π) ./ res .+ (0, 1)))  # CAR has pixels on poles
    end
    nx, ny = shape
    resx, resy = res

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
