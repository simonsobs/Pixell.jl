
"""
    fullsky_geometry([P=CarClenshawCurtis], res; shape = nothing, dims = ())

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
function fullsky_geometry(P::Type{<:CarClenshawCurtis}, res; shape = nothing, dims = ())
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


"""
    pix2sky(m::Enmap, pixcoords)

Convert 1-indexed pixels to sky coordinates. The output sky coordinates are determined by WCS,
but usually are in units of degrees. 

# Arguments:
- `m::Enmap`: the map to obtain the coordinates from
- `pixcoords`: pixcoords should be a 2-d array where "pixcoords[:, i]" is the i-th set of coordinates, 
    or a 1-d array representing a single set of coordinates. 

# Returns: 
- `Array`: same shape as pixcoords

# Examples
```julia-repl
julia> shape, wcs = fullsky_geometry(deg2rad(1))
       m = Enmap(rand(shape...), wcs)
julia> pix_to_world(m, [1.0, 1.0])  # 1-indexing
2-element Vector{Float64}:
 180.0
 -90.0
```
"""
function pix2sky end 

# The default fallback is to call WCS, which calls the C library.
# If you wanted to replicate Pixell behavior, add 1 to x and y of pixcoords.
pix2sky(m::Enmap, pixcoords) = pix_to_world(getwcs(m), pixcoords)
pix2sky!(m::Enmap, pixcoords, skycoords) = pix_to_world!(getwcs(m), pixcoords, skycoords)

# custom Julia-only implementation for one coordinate
function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 pixcoords::Base.AbstractVecOrTuple) where {T,N,AA<:AbstractArray{T,N}}
    
    # retrieve WCS info
    w = getwcs(m)
    α₀, δ₀ = crval(w)
    iα₀, iδ₀ = crpix(w)
    Δα, Δδ = cdelt(w)

    # compute RA (α) and DEC (δ)
    iα, iδ = pixcoords
    α = α₀ + (iα - iα₀) * Δα
    δ = δ₀ + (iδ - iδ₀) * Δδ

    return α, δ
end

function pix2sky!(m::Enmap{T,N,AA,CarClenshawCurtis}, 
                 pixcoords::AbstractArray{Tp,2}, skycoords::AbstractArray{Tp,2}
                 ) where {T,N,AA<:AbstractArray{T,N},Tp}
    
    # retrieve WCS info
    w = getwcs(m)
    α₀, δ₀ = crval(w)
    iα₀, iδ₀ = crpix(w)
    Δα, Δδ = cdelt(w)

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

function pix2sky(m::Enmap{T,N,AA,CarClenshawCurtis}, 
        pixcoords::AbstractArray{Tp,2}) where {T,N,AA<:AbstractArray{T,N},Tp}
    skycoords = similar(pixcoords)
    return pix2sky!(m, pixcoords, skycoords)
end


sky2pix(m::Enmap, skycoords) = world_to_pix(getwcs(m), skycoords)
sky2pix!(m::Enmap, skycoords, pixcoords) = world_to_pix!(getwcs(m), skycoords, pixcoords)
