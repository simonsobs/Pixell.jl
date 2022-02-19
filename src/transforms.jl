    # sharp_make_cc_geom_info_stripe(tot_sphere_rings,
	# 			   cs->nx, cs->phi0,
	# 			   1,cs->nx, //stride_lon ,stride_lat
	# 			   &geom_info,
	# 			   cs->ny,first_ring); //nsubrings, start_index of ring



"""Number of pixels in full CAR ring of this WCS."""
ringsize(wcs::AbstractWCSTransform) = 
    round(Int64, abs(2π / (get_unit(wcs) * cdelt(wcs)[1])))

"""Number of rings in a fullsky version of this WCS."""
ringnum(wcs::AbstractWCSTransform) = 
    1 + round(Int64, abs(π / (get_unit(wcs) * cdelt(wcs)[2])))


# TODO: make this cut everything, and make a full ring expander
function sharp_make_cc_geom_info_stripe(nrings::I, ppring::I, phi0::Real, stride_lon::I, 
                                        stride_lat::I, nsubrings::I, i0::I) where {I<:Integer}
    geom_info_ptr = Ref{Ptr{Cvoid}}()
    μ = chebyshevjacobimoments1(Float64, nrings, 0.0, 0.0)
    weights = Cdouble.(clenshawcurtisweights(μ) .* 2π ./ ppring)
    theta = Cdouble.(range(π, 0.0, length=nrings))
    phi0s = fill(Cdouble(phi0), nrings)
    stride = ones(Cint, nrings)
    offsets = Cptrdiff_t.(ppring .* (0:(nrings-1)))
    nph = fill(Cint(ppring), nrings)

    ccall(
        (:sharp_make_geom_info, Libsharp.libsharp2),
        Cvoid,
        (Cint, Ref{Cint}, Ref{Cptrdiff_t}, Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Ptr{Cvoid}}),
        Cint(nrings),nph, offsets,   stride,    phi0s,         theta,        weights,      geom_info_ptr
    )

    GeomInfo(geom_info_ptr[])
end
