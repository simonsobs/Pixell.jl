
# arbitrary WCS uses the WCSTransform structure from WCS.jl

# retrieve nonallocating WCS information. returns tuples of cdelt, crpix, crval
function unsafe_wcs_read_pair(wcs, sym::Symbol)
    field = WCS.getfield(wcs, sym)
    return unsafe_load(field, 1), unsafe_load(field, 2)
end
getcdelt(wcs::WCSTransform) = unsafe_wcs_read_pair(wcs, :cdelt)
getcrpix(wcs::WCSTransform) = unsafe_wcs_read_pair(wcs, :crpix)
getcrval(wcs::WCSTransform) = unsafe_wcs_read_pair(wcs, :crval)


# internal function for getting the unit of the angles in the WCS header
# multiply the result of this function with your angle to get radians.
getunit(wcs::AbstractWCSTransform) = getunit(Float64, wcs)
function getunit(T::Type{<:Real}, wcs::WCSTransform)
    cunit1, cunit2 = wcs.cunit
    @assert cunit1 == cunit2 "Units of RA and DEC must be the same."
    cunit = cunit1
    if cunit == "deg"
        return T(π/180)
    elseif cunit == "rad"
        return one(T)
    elseif cunit == "arcmin"
        return T(π/180/60)
    elseif cunit == "arcsec"
        return T(π/180/60/60)
    elseif cunit == "mas"
        return T(π/180/60/60/1000)
    end
    @warn "Can't recognize the WCS unit: $(cunit). Assuming degrees."
    return T(π/180)  # give up and return degrees
end

# Select the first n axes, should move to WCS.jl at some point
function sub(src::WCS.WCSTransform, n::Int)
    dst = WCSTransform(n)
    nsub = Array{Cint}([n])
    axes = Array{Cint}(collect(1:n))
    ccall((:wcssub, libwcs), Cint,
          (Cint, Ref{WCSTransform}, Ptr{Cint}, Ptr{Cint}, Ref{WCSTransform}),
          0, src, nsub, axes, dst)
    dst
end

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


"""Check if a WCS is a cylindrical pixelization"""
function iscyl(wcs::WCSTransform)
    # right now we only support CAR
    ctype = wcs.ctype
    if ctype[1] == "RA---CAR" && ctype[2] == "DEC--CAR"
        return true
    end
    return false
end

"""Area of a patch given shape and WCS, in steradians."""
function skyarea(shape, wcs::WCSTransform)
    if iscyl(wcs)
        return skyarea_cyl(shape, wcs)
    end 
    error("Area for this WCS is not implemented.")
end

# generic cylindrical area calculation (user is promising it's cylindrical)
function skyarea_cyl(shape, wcs::AbstractWCSTransform)
    N_α, N_δ = shape
    α, δ = pix2sky(shape, wcs, SVector(0., 0.), SVector(0.5, N_δ+0.5); safe=false)
    δ₁, δ₂ = sort(δ)
    δ₁, δ₂ = max(-π/2, δ₁), min(π/2, δ₂)
    Δα, Δδ = getcdelt(wcs) .* getunit(wcs)
    return (sin(δ₂) - sin(δ₁)) * abs(Δα) * N_α
end

function sliced_wcs(wcs::WCSTransform, cdelt′, crpix′)
    new_wcs = deepcopy(wcs)
    new_wcs.cdelt = collect(cdelt′)
    new_wcs.crpix = collect(crpix′)
    return new_wcs
end
