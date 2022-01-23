using Base: ViewIndex, @propagate_inbounds, AbstractCartesianIndex

abstract type MapProjection end
abstract type EquiCylProjection <: MapProjection end  # equidistant cylindrical projection

struct CarProjection <: EquiCylProjection end  # plate carrée
struct CeaProjection <: EquiCylProjection end  # cylindrical equal area

"""
Map type, contains an AbstractArray and a WCS object, but behaves like the
AbstractArray it contains for array operations.
It only implements the subset of Base.Array operations which are common on maps.
You should work with the data directly using `enmap_instance.data` if you need
additional Array functions.
"""
struct Enmap{T,N,AA<:AbstractArray,P<:MapProjection} <: AbstractArray{T,N}
    data::AA  # some kind of abstract array
    wcs::WCSTransform  # WCS object from WCS.jl
end


function Enmap(data::A, wcs, ::Type{PT}) where {A<:AbstractArray,PT}
    Enmap{eltype(A),ndims(A),A,PT}(data, wcs)
end
# create CAR maps by default
function Enmap(data::A, wcs) where {A<:AbstractArray}
    Enmap{eltype(A),ndims(A),A,CarProjection}(data, wcs)
end

Base.parent(x::Enmap) = x.data
getwcs(x::Enmap) = x.wcs

# forward all array traits to parent. based on MetaArrays
Base.size(x::Enmap) = size(parent(x))
Base.axes(x::Enmap) = Base.axes(parent(x))
Base.IndexStyle(x::Enmap) = IndexStyle(parent(x))

@propagate_inbounds Base.view(A::Enmap, idxs::ViewIndex...) = Enmap(view(parent(A), idxs...), getwcs(A))
@propagate_inbounds Base.view(A::Enmap, idxs::Union{ViewIndex,AbstractCartesianIndex}...) = Enmap(view(parent(A), idxs...), getwcs(A))
@propagate_inbounds Base.view(A::Enmap, idxs...) = Enmap(view(parent(A), idxs...), getwcs(A))

@propagate_inbounds Base.getindex(x::Enmap, i::Int...) = getindex(parent(x), i...)
@propagate_inbounds Base.getindex(x::Enmap, i...) = enmapwrap(x, getindex(parent(x), i...))
@propagate_inbounds Base.setindex!(x::Enmap, v, i...) = (setindex!(parent(x), v, i...); x)

function Base.similar(x::Enmap, ::Type{S}, dims::NTuple{<:Any,Int}) where {S}
    Enmap(getwcs(x), similar(parent(x), S, dims))
end

enmapwrap(x::Enmap{T}, val::T) where {T} = val
enmapwrap(x::Enmap{T,N,AA,P}, val::AbstractArray{T,N}) where {T,N,AA,P} = Enmap{T,N,AA,P}(val, getwcs(x))

# if array slicing ends up changing dimension, follow the parent
function enmapwrap(x::Enmap{T,N,AA,P}, val::AAV) where {T,N,AA,P,NV,AAV<:AbstractArray{T,NV}}
    Enmap{T,NV,AAV,P}(val, getwcs(x))
end
# enmapwrap(x::Enmap, val) = error("Unexpected result type $(typeof(val)).")

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


struct NoWCS end
combine(x::NoWCS, y) = y
combine(x, y::NoWCS) = x
combine(x::NoWCS, ::NoWCS) = x
combine(x::Enmap, y::Enmap) = x  # TODO: check compatibility


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
    Enmap(copy(bc_), bc.args[1][1])
end
