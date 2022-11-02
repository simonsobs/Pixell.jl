"""
Fast custom WCS structure.
"""
struct Gnomonic{T} <: AbstractWCSTransform
    cdelt::Tuple{T,T}
    crpix::Tuple{T,T}
    crval::Tuple{T,T}
    unit::T  # conversion factor to radians
end
getunit(T::Type{<:Real}, wcs::Gnomonic) = T(wcs.unit)
getcdelt(wcs::Gnomonic) = wcs.cdelt
getcrpix(wcs::Gnomonic) = wcs.crpix
getcrval(wcs::Gnomonic) = wcs.crval

Base.copy(w::W) where {W<:Gnomonic} = w  # Gnomonic is fully immutable

function Base.convert(::Type{Gnomonic{T}}, w0::WCSTransform) where T
    return Gnomonic{T}(
        T.(getcdelt(w0)), T.(getcrpix(w0)), T.(getcrval(w0)), getunit(T, w0))
end

function Base.convert(::Type{WCSTransform}, w0::Gnomonic{T}) where T
    return WCSTransform(2;
        ctype = ["RA---TAN", "DEC--TAN"],
        cdelt = collect(getcdelt(w0)),
        crpix = collect(getcrpix(w0)),
        crval = collect(getcrval(w0)))
end

# this kind of WCS only has two spatial dimensions. this check should be constant-propagated
function Base.getproperty(wcs::Gnomonic, k::Symbol)
    if k == :naxis
        return 2
    end
    return getfield(wcs, k)
end

function Base.show(io::IO, wcs::Gnomonic{T}) where T
    expr = join(["$k=$(getproperty(wcs, Symbol(k)))"
                 for k in ["naxis","cdelt","crval","crpix"]], ",")
    print(io, "Gnomonic{$(T)}($expr)")
end

function sky2pix(shape, wcs::Gnomonic{T}, ra, dec; safe=false) where T
    scale = one(T) / first(wcs.cdelt)
    unit = getunit(wcs)
    α₀, δ₀ = deg2rad.(wcs.crval)
    α, δ = ra, dec
    A = cos(δ) * cos(α - α₀)
    F = scale / unit / (sin(δ₀) * sin(δ) + A * cos(δ₀))
    
    LINE = -F * ( cos(δ₀) * sin(δ) - A * sin(δ₀) )
    SAMPLE = -F * cos(δ) * sin(α - α₀)
    
    X, Y = wcs.crpix .- (SAMPLE, LINE)
    return X, Y
end

function pix2sky(shape, wcs::Gnomonic{T}, ra_pixel, dec_pixel; safe=false) where T
    scale = one(T) / first(wcs.cdelt)
    unit = getunit(wcs)
    α₀, δ₀ = deg2rad.(wcs.crval)
    X = (ra_pixel - wcs.crpix[1]) * unit / scale
    Y = (dec_pixel - wcs.crpix[2]) * unit / scale
    
    D = atan(√(X^2 + Y^2))
    B = atan(-X/Y)
    XX = sin(δ₀) * sin(D) * cos(B) + cos(δ₀) * cos(D)
    YY = sin(D) * sin(B)
    
    α = α₀ + atan(YY/XX)
    δ = asin(sin(δ₀) * cos(D) - cos(δ₀) * sin(D) * cos(B))

    return α, δ
end
