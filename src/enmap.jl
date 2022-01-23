
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

Base.parent(x::Enmap) = x.data

# forward all array traits to parent
Base.size(x::Enmap) = size(parent(x))
Base.axes(x::Enmap) = Base.axes(parent(x))
Base.IndexStyle(x::Enmap) = IndexStyle(parent(x))

Base.@propagate_inbounds Base.getindex(x::Enmap, i::Int...) = getindex(parent(x), i...)
