

using Pixell
using Plots
using StaticArrays

##

function distance_transform_bruteforce(x)
    ax1 = 2:(size(x,1)-1)
    ax2 = 2:(size(x,2)-1)
    bads = Tuple{Int,Int}[]
    for j in ax2, i in ax1
        if iszero(x[i,j])
            push!(bads, (i,j))
        end
    end
    result = fill(0.0, size(x))
    for j in ax2, i in ax1
        min_dist = MAX_DIST
        for (bi, bj) in bads
            this_dist = metric(bi-i, bj-j, 0, 0)
            min_dist = min(min_dist, this_dist)
        end
        result[i,j] = min_dist
    end
    return sqrt.(result)
end




##

##

# vecs = fill(Int(99), (2, size(x)...) )


function metric(x1, y1, x2, y2)
    return (x1 - x2)^2 + (y1 - y2)^2
    # i′ = i+iof
    # j′ = j+jof
    # return (v[1,i′,j′] - iof)^2 + (v[2,i′,j′] - jof)^2
end

const mask_1a = ((-1,-1), (0, -1), (1,-1), (-1,0), (0,0))
const mask_1b = ((0,0), (1,0))
const mask_2a = ((-1,0), (0,0))
const mask_2b = ((0,0), (1,0), (-1,1), (0, 1), (1,1))


##

##

function op!(vecs, i, j, mask)

    min_dist = MAX_DIST

    # loop over all vectors at positions in the mask
    for (i_mask, (iof, jof)) in enumerate(mask)
        i′ = i+iof
        j′ = j+jof
        for v in vecs[i′,j′]
            this_dist = metric(v[1] - iof, v[2] - jof, 0, 0)
            if this_dist < min_dist
                min_dist = this_dist
            end
        end
    end

    # can check 0,0 here
    if min_dist > 0
        newvecs = SVector{2, Int}[]
        for (iof, jof) in mask
            for v in vecs[i+iof,j+jof]
                this_dist = metric(v[1] - iof, v[2] - jof, 0, 0)
                if this_dist <= min_dist #+ 1/2
                    xv = SVector{2, Int}(v[1] - iof, v[2] - jof)
                    if xv ∉ newvecs
                        push!(newvecs, xv)
                    end
                end
            end
        end
        vecs[i,j] = newvecs
    end
    # @show newvecs
end


function do_ops(vecs)
    for j in ax2, i in ax1
        op!(vecs, i, j, mask_1a)
    end

    for j in ax2, i in reverse(ax1)
        op!(vecs, i, j, mask_1b)
    end

    for j in reverse(ax2), i in ax1
        op!(vecs, i, j, mask_2a)
    end

    for j in reverse(ax2), i in reverse(ax1)
        op!(vecs, i, j, mask_2b)
    end
end


##


x = ones(128,128)


const MAX_DIST = 10*prod(size(x))
const ε = 4.0


for i in 1:30
    x[rand(2:(size(x,1)-1)), rand(2:(size(x,2)-1))] = 0.0
end

vecs = Matrix{Vector{SVector{2, Int}}}(undef, size(x))

for j in axes(vecs,2), i in axes(vecs,1)
    if iszero(x[i,j])
        vecs[i,j] = SVector{2, Int}[SA[0, 0]]
    else 
        vecs[i,j] = SVector{2, Int}[SA[MAX_DIST, MAX_DIST]]
    end
end

# x[12,6] = 0.0
# x[5,18] = 0.0
# x[8,12] = 0.0
ax1 = 2:(size(x,1)-1)
ax2 = 2:(size(x,2)-1)
heatmap(x', aspectratio=1)

@time do_ops(vecs)

distmap = zeros(size(x))
for j in ax2, i in ax1
    min_dist = MAX_DIST
    for v in vecs[i,j]
        this_dist = metric(v[1], v[2], 0, 0)
        min_dist = min(min_dist, this_dist)
    end
    distmap[i,j] = sqrt(min_dist)
end

heatmap(distmap', aspectratio=1, clim=(0,sqrt(prod(size(x)))/2))


##
distmap_bf = distance_transform_bruteforce(x)

##
heatmap(distmap_bf .- distmap, aspectratio=1)

##

sum(abs.(distmap_bf .- distmap))
