using Base: ViewIndex, @propagate_inbounds, AbstractCartesianIndex

"""
Map type, contains an AbstractArray and a WCS object, but behaves like the
AbstractArray it contains for array operations.
It only implements the subset of Base.Array operations which are common on maps.
You should work with the data directly using `enmap_instance.data` if you need
additional Array functions.
"""
struct Enmap{T,N,AA<:AbstractArray,W<:AbstractWCSTransform} <: AbstractArray{T,N}
    data::AA  # some kind of abstract array
    wcs::W    # WCS object
end


function Enmap(data::A, wcs::W) where {A<:AbstractArray,W<:AbstractWCSTransform}
    Enmap{eltype(A),ndims(A),A,W}(data, wcs)
end

struct NoWCS end

Base.parent(x::Enmap) = x.data
getwcs(x::Enmap) = x.wcs
getwcs(x) = NoWCS()


function Base.show(io::IO, imap::Enmap)
    expr = "Enmap(shape=$(size(imap)),wcs=$(imap.wcs))"
    print(io, expr)
end


# forward all array traits to parent. based on MetaArrays
Base.size(x::Enmap) = size(parent(x))
Base.axes(x::Enmap) = Base.axes(parent(x))
Base.IndexStyle(x::Enmap) = IndexStyle(parent(x))

const NotIntIndex = Union{AbstractArray,AbstractCartesianIndex,Colon}
const AnyIndex = Union{ViewIndex,AbstractCartesianIndex,Colon}
@propagate_inbounds function Base.view(x::Enmap, idxs...) 
    new_shape, new_wcs = slice_geometry(x, idxs...)
    Enmap(view(parent(x), idxs...), new_wcs)
end


# totally give up if the RA and DEC slices are eliminated, and just return a view of the parent
@propagate_inbounds Base.view(x::Enmap, ix::Int, iy::NotIntIndex, idxs::AnyIndex...) = 
    view(parent(x), ix, iy, idxs...)
@propagate_inbounds Base.view(x::Enmap, ix::NotIntIndex, iy::Int, idxs::AnyIndex...) = 
    view(parent(x), ix, iy, idxs...)
@propagate_inbounds Base.view(x::Enmap, ix::Int, iy::Int, idxs::AnyIndex...) = 
    view(parent(x), ix, iy, idxs...)

@propagate_inbounds Base.getindex(x::Enmap, i::Int...) = getindex(parent(x), i...)
@propagate_inbounds function Base.getindex(x::Enmap, i...)
    enmapwrap(x, getindex(parent(x), i...), i...)
end
@propagate_inbounds Base.setindex!(x::Enmap, v, i...) = (setindex!(parent(x), v, i...); x)

function Base.similar(x::Enmap, ::Type{S}, dims::NTuple{<:Any,Int}) where {S}
    Enmap(similar(parent(x), S, dims), deepcopy(getwcs(x)))
end

enmapwrap(x::Enmap{T}, val::T, i...) where {T} = val
function enmapwrap(x::Enmap{T,N,AA,W}, val::AbstractArray{T,N}, i...) where {T,N,AA,W} 
    new_shape, new_wcs = slice_geometry(x, i...)
    Enmap{T,N,AA,W}(val, new_wcs)
end

# if array slicing ends up changing dimension, follow the parent
function enmapwrap(x::Enmap{T,N,AA,W}, val::AAV, i...) where {T,N,AA,W,NV,AAV<:AbstractArray{T,NV}}
    new_shape, new_wcs = slice_geometry(x, i...)
    Enmap{T,NV,AAV,W}(val, new_wcs)
end

enmapwrap(x::Enmap{T,N,AA,W}, val::AAV, i1::Int, i...) where {T,N,AA,W,NV,AAV<:AbstractArray{T,NV}} = val
enmapwrap(x::Enmap{T,N,AA,W}, val::AAV, i1, i2::Int, i...) where {T,N,AA,W,NV,AAV<:AbstractArray{T,NV}} = val
enmapwrap(x, val, i...) = error("Unexpected result type $(typeof(val)).")

Base.strides(x::Enmap) = strides(parent(x))
Base.stride(x::Enmap, i::Int) = stride(parent(x), i)


# the meta array broadcast style should retain the nested style information for
# whatever array type the meta array wraps
struct EnmapStyle{S} <: Broadcast.BroadcastStyle end
EnmapStyle(s::S) where {S<:Broadcast.BroadcastStyle} = EnmapStyle{S}()
Base.Broadcast.BroadcastStyle(::Type{<:Enmap{T,N,A}}) where {T,N,A} =
    enmapstyle(Broadcast.BroadcastStyle(A))
Base.Broadcast.BroadcastStyle(a::EnmapStyle{A}, b::EnmapStyle{B}) where {A,B} =
    enmapstyle(Broadcast.BroadcastStyle(A(), B()))
function Base.Broadcast.BroadcastStyle(a::EnmapStyle{A}, b::B) where
{A,B<:Broadcast.BroadcastStyle}

    a_ = A()
    left = enmapstyle(Broadcast.BroadcastStyle(a_, b))
    if !(left isa Broadcast.Unknown)
        left
    else
        enmapstyle(Broadcast.BroadcastStyle(b, a_))
    end
end
enmapstyle(x) = EnmapStyle(x)
enmapstyle(x::Broadcast.Unknown) = x


combine(x::NoWCS, y) = copy(y)
combine(x, y::NoWCS) = copy(x)
combine(x::NoWCS, ::NoWCS) = x
combine(x::AbstractWCSTransform, y::AbstractWCSTransform) = copy(x)  # TODO: check compatibility


####
#### broadcasting wrap from MetaArrays
function enmap_broadcasted(bc::Broadcast.Broadcasted{S}, wcss) where {S}
    args = enmap_.(bc.args, wcss)
    Broadcast.Broadcasted{EnmapStyle{S}}(bc.f, args, bc.axes)
end
enmap_broadcasted(result, wcss) = Enmap(result, reduce(combine, wcss))

enmap_(x, ::NoWCS) = x
enmap_(x, wcs) = Enmap(x, wcs)
maybeparent_(x) = x
maybeparent_(x::Enmap) = parent(x)
getwcs_(x) = NoWCS()
getwcs_(x::Enmap) = getwcs(x)

# broadcasted:
function Base.Broadcast.broadcasted(::EnmapStyle{S}, f, xs...) where {S}
    bc = Broadcast.broadcasted(S(), f, maybeparent_.(xs)...)
    enmap_broadcasted(bc, getwcs_.(xs))
end

# instantiate:
# after instantiation, the broadcasted object is flattened and the
# argument contains all meteadata
function Base.Broadcast.instantiate(bc::Broadcast.Broadcasted{M}) where {S,M<:EnmapStyle{S}}

    # simplify
    bc_ = Broadcast.flatten(bc)
    # instantiate the nested broadcast (that the meta array wraps)
    bc_nested = Broadcast.Broadcasted{S}(bc_.f, maybeparent_.(bc_.args))
    inst = Broadcast.instantiate(bc_nested)
    # extract and combine the meta data
    meta = reduce(combine, getwcs_.(bc_.args))
    # place the meta data on the first argument
    args = ((meta, inst.args[1]), Base.tail(inst.args)...)
    # return the instantiated metadata broadcasting
    Broadcast.Broadcasted{M}(bc_.f, args, bc_.axes)
end

# similar: becuase we bypass the default broadcast machinery no call to similar
# is made directly: similar will be called within the nested call to `copy`
# after the metadata has been stripped and the wrapped array is exposed

# copyto!:
function Base.copyto!(dest::AbstractArray, bc::Broadcast.Broadcasted{<:EnmapStyle{S}}) where {S}
    args_ = (bc.args[1][2], Base.tail(bc.args)...)
    bc_ = Broadcast.Broadcasted{S}(bc.f, args_, bc.axes)
    copyto!(dest, bc_)
end

function Base.copyto!(dest::Enmap, bc::Broadcast.Broadcasted{Nothing})
    copyto!(parent(dest), bc)
end

# copy:
function Base.copy(bc::Broadcast.Broadcasted{<:EnmapStyle{S}}) where {S}
    # because the axes have been instantiated, we can safely assume the first
    # argument contains the meta data
    args_ = (bc.args[1][2], Base.tail(bc.args)...)
    bc_ = Broadcast.Broadcasted{S}(bc.f, args_, bc.axes)
    Enmap(copy(bc_), copy(bc.args[1][1]))
end

Base.deepcopy(x::Enmap) = Enmap(deepcopy(parent(x)), deepcopy(getwcs(x)))

function resolve_polcconv!(data::A, header::FITSIO.FITSHeader, sel; verbose=true) where {A<:AbstractArray}
    naxis = header["NAXIS"]
    has_sel = length(sel) == naxis
    for i in 1:naxis
        if getkey(header, "CTYPE$i", "") == "STOKES"
            signs_size = ones(Int, header["NAXIS"])
            signs_size[i] = 3
            verbose && (println("convert to IAU: flip U in axis $i"))
            signs = ones(eltype(A), signs_size...)
            signs[signs_size] .= -1
            if has_sel
                data .*= selectdim(signs, i, sel[i])
            else
                data .*= signs
            end
        end
    end
end

# read fits file into Enmap
function read_map(path::String; hdu::Int=1, sel=(), wcs::Union{WCSTransform,Nothing}=nothing, verbose=true, trim=true)
    f = FITS(path, "r")
    data = read(f[hdu], sel...)
    if isnothing(wcs)
        # parse header from hdu as FITSIO.FITSHeader
        header = read_header(f[hdu])
        # handle IAU <--> COSMOS conversion
        if "STOKES" in header.values
            # default to no flipping (COSMO)
            verbose && !haskey(header, "POLCCONV") && (println("STOKES found but POLCCONV not found, assuming IAU"))
            polcconv = getkey(header, "POLCCONV", "COSMO")
            polcconv == "IAU" && (resolve_polcconv!(data, header, sel; verbose=verbose))
        end
        # WCS.from_header expects each key to be right-padded with len=80
        header_str = join([@sprintf("%-80s", f) for f in split(string(header), "\n")])
        if trim == false
            wcs = sub(WCS.from_header(header_str)[1], 2)
        else
            wcs0 = WCS.from_header(header_str)[1]
            @assert wcs0.ctype[1] == "RA---CAR" 
            @assert wcs0.ctype[2] == "DEC--CAR"
            wcs = convert(CarClenshawCurtis{Float64}, wcs0)
        end
    end
    Enmap(data, wcs)
end

function write_map(fname::String, emap::Enmap)
    f = FITS(fname, "w")
    # it is important to write data first. The data will be written with
    # some default headers that describe the shape for us. We then update
    # header based on wcs information to populate the rest of the cards
    write(f, emap.data)
    header = WCS.to_header(Base.convert(WCSTransform, getwcs(emap)))
    cards = [header[1+(i-1)*80:i*80] for i = 1:round(Int, length(header)/80)]
    for c in cards
        FITSIO.CFITSIO.fits_write_record(f.fitsfile, c)
    end
    close(f)
end
