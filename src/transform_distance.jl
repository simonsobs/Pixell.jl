

"""Distance Transform (DT), especially Spherical Distance Transforms (SDT)"""
abstract type AbstractDT end
abstract type AbstractSDT <: AbstractDT end
abstract type AbstractSeqSDT <: AbstractSDT end
struct BruteForceSDT <: AbstractSDT end
struct ApproxSeqSDT <: AbstractSeqSDT end

struct ExactSeqSDT <: AbstractSeqSDT
    ϵ::Float64
    buffer::Vector{Tuple{Int,Int}}
end
ExactSeqSDT(ϵfactor=1.0) = ExactSeqSDT(ϵfactor, Tuple{Int,Int}[])


const DANIELSSON_MASK1A = ((-1,-1), (0, -1), (1,-1), (-1,0), (0,0))
const DANIELSSON_MASK1B = ((0,0), (1,0))
const DANIELSSON_MASK2A = ((0,0), (1,0), (-1,1), (0, 1), (1,1))
const DANIELSSON_MASK2B = ((-1,0), (0,0))


"""Brute-force spherical distance transform. ~ O(Nₚᵢₓ × N₀). Doesn't handle wrapping."""
function distance_transform(DT::BruteForceSDT, m::Enmap)
    αs = pix2sky(m, collect(1:size(m,1)), ones(size(m,1)))[1]
    δs = pix2sky(m, ones(size(m,2)), collect(1:size(m,2)))[2]
    
    ax1 = 2:(size(m,1)-1)
    ax2 = 2:(size(m,2)-1)


    bads = Tuple{Int,Int}[]
    for j in ax2, i in ax1
        if iszero(m[i,j])
            push!(bads, (i,j))
        end
    end
    result = fill(0.0, size(m))
    
    max_dist = Inf
    for j in ax2, i in ax1
        min_dist = max_dist
        for (bi, bj) in bads
            this_dist = metric_from_indices(DT, αs, δs, i, j, bi, bj)
            min_dist = min(min_dist, this_dist)
        end
        result[i,j] = min_dist
    end
    return Enmap(sqrt.(result), getwcs(m))
end

"""Metric on the sphere for RA (α) and DEC (δ)"""
function metric(::AbstractSDT, α₁, δ₁, α₂, δ₂)
    x₁ = cos(δ₁) * cos(α₁)
    y₁ = cos(δ₁) * sin(α₁) 
    z₁ = sin(δ₁)

    x₂ = cos(δ₂) * cos(α₂)
    y₂ = cos(δ₂) * sin(α₂) 
    z₂ = sin(δ₂)
    res = (x₁ - x₂)^2 + (y₁ - y₂)^2 + (z₁ - z₂)^2
    return res
end

# we initialize the vectors as something huge, and return a huge metric for that initial condition
function metric_from_indices(DT::AbstractSDT, αs, δs, i₁, j₁, i₂, j₂)
    if (checkbounds(Bool, αs, i₁) && checkbounds(Bool, δs, j₁) 
            && checkbounds(Bool, αs, i₂) && checkbounds(Bool, δs, j₂))
        return metric(DT, αs[i₁], δs[j₁], αs[i₂], δs[j₂])
    end
    return Inf
end


function propagate!(DT::ApproxSeqSDT, αs, δs, v, i, j, mask)

    min_dist = Inf
    i_min = 0
    j_min = 0
    for (iof, jof) in mask
        i′ = i + iof
        j′ = j + jof
        dist = metric_from_indices(DT, αs, δs, 
            i + v[1,i′,j′] + iof, j + v[2,i′,j′] + jof, 
            i, j)
        if dist < min_dist
            min_dist = dist 
            i_min =  v[1,i′,j′] + iof
            j_min =  v[2,i′,j′] + jof
        end
    end
    v[1,i,j] = i_min
    v[2,i,j] = j_min
end


function setup_distance_vectors(::ApproxSeqSDT, m::Enmap)
    shape = size(m)
    max_vector_comp = typemax(Int)
    v = fill(max_vector_comp, (2, shape...) )

    ax1 = 2:(size(m,1)-1)
    ax2 = 2:(size(m,2)-1)
    for j in ax2, i in ax1
        if iszero(m[i,j])
            v[:,i,j] .= (0, 0)
        end
    end
    return v
end

function distance_transform_vectors(DT::AbstractSeqSDT, m::Enmap)
    αs = pix2sky(m, collect(1:size(m,1)), ones(size(m,1)))[1]
    δs = pix2sky(m, ones(size(m,2)), collect(1:size(m,2)))[2]
    v = setup_distance_vectors(DT, m)

    ax1 = 2:(size(m,1)-1)
    ax2 = 2:(size(m,2)-1)

    for j in ax2
        for i in ax1
            propagate!(DT, αs, δs, v, i, j, DANIELSSON_MASK1A)
        end
        for i in reverse(ax1)
            propagate!(DT, αs, δs, v, i, j, DANIELSSON_MASK1B)
        end
    end

    for j in reverse(ax2)
        for i in reverse(ax1)
            propagate!(DT, αs, δs, v, i, j, DANIELSSON_MASK2A)
        end
        for i in (ax1)
            propagate!(DT, αs, δs, v, i, j, DANIELSSON_MASK2B)
        end
    end

    return v, αs, δs
end

function distance_transform(DT::ApproxSeqSDT, m::Enmap)
    # these vectors point from the pixel to the offending pixel that it's closest to
    v, αs, δs = distance_transform_vectors(DT::ApproxSeqSDT, m::Enmap)
    distmap = zeros(size(m))
    for j in 2:(size(distmap,2)-1), i in 2:(size(distmap,1)-1)
        i′, j′ = i + v[1,i,j], j + v[2,i,j]
        distmap[i,j] = sqrt.(metric_from_indices(DT, αs, δs, i, j, i′ ,j′))
    end
    return Enmap(distmap, Pixell.getwcs(m))
end

function setup_distance_vectors(::ExactSeqSDT, m::Enmap)
    shape = size(m)
    max_vector_comp = typemax(Int)
    v = Array{Vector{Tuple{Int,Int}},2}(undef, shape)

    for j in axes(m, 2), i in axes(m, 1)
        v[i,j] = Vector{Tuple{Int,Int}}[]
        if iszero(m[i,j])
            push!(v[i,j], (0,0))
        else
            push!(v[i,j], (max_vector_comp, max_vector_comp))
        end
    end
    return v
end

function propagate!(DT::ExactSeqSDT, αs, δs, vecs, i, j, mask)

    min_dist = Inf
    ϵ = metric_from_indices(DT, αs, δs, 1, 1, 2, 2) * DT.ϵ
    # loop over all vectors at positions in the mask
    for (iof, jof) in mask
        i′ = i+iof
        j′ = j+jof
        for (v1, v2) in vecs[i′,j′]
            this_dist = metric_from_indices(DT, αs, δs, 
                i + v1 + iof, j + v2 + jof, 
                i, j)
            if this_dist < min_dist
                min_dist = this_dist
            end
        end
    end
    
    if isfinite(min_dist) && min_dist > 0
        empty!(DT.buffer)
        for (iof, jof) in mask
            for (v1, v2) in vecs[i+iof,j+jof]
                this_dist = metric_from_indices(DT, αs, δs, 
                    i + v1 + iof, j + v2 + jof, 
                    i, j)
                if this_dist <= min_dist + ϵ
                    xv = (v1 + iof, v2 + jof)
                    if xv ∉ DT.buffer
                        push!(DT.buffer, xv)
                    end
                end
            end
        end
        resize!(vecs[i,j], length(DT.buffer))
        vecs[i,j] .= DT.buffer
    end
end


function distance_transform(DT::ExactSeqSDT, m::Enmap)
    # these vectors point from the pixel to the offending pixel that it's closest to
    v, αs, δs = distance_transform_vectors(DT, m)
    distmap = zeros(size(m))
    max_dist = Inf

    for j in 2:(size(distmap,2)-1), i in 2:(size(distmap,1)-1)
        
        min_dist = max_dist
        for v in v[i,j]
            i′, j′ = i + v[1], j + v[2]
            this_dist = sqrt.(metric_from_indices(DT, αs, δs, i, j, i′ ,j′))
            min_dist = min(min_dist, this_dist)
        end
        distmap[i,j] = min_dist
    end
    return Enmap(distmap, Pixell.getwcs(m))
end
##

##