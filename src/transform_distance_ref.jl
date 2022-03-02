

using Pixell
using Plots

x = ones(64,64)
x[13,48] = 0.0
x[22,22] = 0.0
x[44,62] = 0.0
heatmap(x', aspectratio=1)

##
const MAX_DIST = prod(size(x))
const ax1 = 2:(size(x,1)-1)
const ax2 = 2:(size(x,2)-1)

vecs = fill(Int(MAX_DIST), (2, size(x)...) )

for j in ax2, i in ax1
    if iszero(x[i,j])
        vecs[:,i,j] .= (0, 0)
    end
end

function dist(v, i, j, iof, jof)
    i′ = i+iof
    j′ = j+jof
    return (v[1,i′,j′] - iof)^2 + 1.0 * (v[2,i′,j′] - jof)^2
end


function op!(v, i, j, mask)
    ind_maskmin = argmin((dist(v, i, j, of[1], of[2]) for of in mask))
    iof, jof = mask[ind_maskmin]
    i′ = i+iof
    j′ = j+jof
    v[:,i,j] .= v[1,i′,j′] - iof, v[2,i′,j′] - jof
end

function do_op(vecs)
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
#     kmin = argmin((dist(vecs, i, j, of[1], of[2]) for of in mask_2b))
#     vecs[:,i,j] .= vecs[:,i,j] .+ mask_2b[kmin]
# end

##

@time do_op(vecs)
# vecs[:,128,48] .= 1000.0
# vecs[:,50,180] .= 1000.0
# vecs[:,180,112] .= 1000.0

distmap = zeros(size(x))
for j in ax2, i in ax1
    distmap[i,j] = sqrt.(vecs[1,i,j].^2 .+ vecs[2,i,j].^2)
end
heatmap(distmap', aspectratio=1)


##
@btime t[2][(t[1],1)] setup=(t=(rand(0:1), 
    Dict{Tuple{Int,Int},Float64}((1,1)=>0.0, (0,1)=>43.0)))

##

