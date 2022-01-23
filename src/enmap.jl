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


function Enmap(data::A, wcs, ::Type{PT}) where {A<:AbstractArray, PT}
    Enmap{eltype(A), ndims(A), A, PT}(data,wcs)
end
# create CAR maps by default
function Enmap(data::A, wcs) where {A<:AbstractArray}
    Enmap{eltype(A),ndims(A),A,CarProjection}(data,wcs)
end

Base.parent(x::Enmap) = x.data
getwcs(x::Enmap) = x.wcs

# forward all array traits to parent. based on MetaArrays
Base.size(x::Enmap) = size(parent(x))
Base.axes(x::Enmap) = Base.axes(parent(x))
Base.IndexStyle(x::Enmap) = IndexStyle(parent(x))

@propagate_inbounds Base.view(A::Enmap, idxs::ViewIndex...) = Enmap(view(parent(A), idxs...), getwcs(A))
@propagate_inbounds Base.view(A::Enmap, idxs::Union{ViewIndex,AbstractCartesianIndex}...) = Enmap(view(parent(A),idxs...), getwcs(A))
@propagate_inbounds Base.view(A::Enmap, idxs...) = Enmap(view(parent(A), idxs...), getwcs(A))

@propagate_inbounds Base.getindex(x::Enmap, i::Int...) = getindex(parent(x), i...)
@propagate_inbounds Base.getindex(x::Enmap, i...) = enmapwrap(x, getindex(parent(x), i...))
@propagate_inbounds Base.setindex!(x::Enmap, v, i...) = (setindex!(parent(x), v, i...); x)

function Base.similar(x::Enmap,::Type{S},dims::NTuple{<:Any,Int}) where S
    Enmap(getwcs(x), similar(parent(x), S, dims))
end

enmapwrap(x::Enmap{T,N,AA,P},val::T) where {T,N,AA,P} = val
enmapwrap(x::Enmap{T,N,AA,P},val::AbstractArray{T,N}) where {T,N,AA,P} = Enmap{T,N,AA,P}(val, getwcs(x))

# if array slicing ends up changing dimension, follow the parent
function enmapwrap(x::Enmap{T,N,AA,P},val::AAV) where {T,N,AA,P,NV,AAV<:AbstractArray{T,NV}} 
    Enmap{T,NV,AAV,P}(val, getwcs(x))
end
# enmapwrap(x::Enmap, val) = error("Unexpected result type $(typeof(val)).")

Base.strides(x::Enmap) = strides(parent(x))
Base.stride(x::Enmap, i::Int) = stride(parent(x), i)
