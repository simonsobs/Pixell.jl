
"""
Map type, contains an AbstractArray and a WCS object, but behaves like the
AbstractArray it contains for array operations.
It only implements the subset of Base.Array operations which are common on maps.
You should work with the data directly using `enmap_instance.data` if you need
additional Array functions.
"""
struct Enmap{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    data::AA  # some kind of abstract array
    wcs::WCSTransform  # WCS object from WCS.jl
end

