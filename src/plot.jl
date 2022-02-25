module Enplot

using Colors: RGB, weighted_color_mean
using RecipesBase
using ColorSchemes
import Pixell: Enmap, pix2sky

# convert 0:255 integer to 0:1 Float64
rgb(c...) = RGB{Float64}((c./255)...)

# build a ColorScheme based on a list of colors and their locations in
# 0:1 colormap
function build_colorscheme(colors, locs)
    xs = collect(0:0.01:1)
    map(xs) do x
        i = min(findlast((<=)(x), locs), length(locs)-1)
        r = (x-locs[i])/(locs[i+1]-locs[i])
        weighted_color_mean(1-r, colors[i], colors[i+1])
    end |> ColorScheme
end

# here we define some common color schemes
cschemes = Dict{Symbol,ColorScheme}()
cschemes[:planck] = build_colorscheme([rgb(0, 0, 255), rgb(0, 215, 255), rgb(255, 237, 217), rgb(255, 180, 0), rgb(255, 75, 0), rgb(100, 0, 0)], [0, 0.332, 0.5, 0.664, 0.828, 1])

function register_colors!()
    for (k, v) in cschemes
        ColorSchemes.colorschemes[k] = v
    end
end

export register_colors!

# here we define common plot recipes
@recipe function f(imap::Enmap)
    seriestype   := :heatmap
    aspect_ratio := :equal
    xformatter   := x -> pix2sky_formatter(x, imap)
    yformatter   := x -> pix2sky_formatter(x, imap; ind=2)
    color          --> :planck
    xlim           --> (1, size(imap.data,1))  # not nice
    ylim           --> (1, size(imap.data,2))
    colorbar       --> true
    minorticks     --> true
    tick_direction --> :out
    grid           --> false
    size           --> (500, 
                        500 * size(imap.data,2) / size(imap.data,1))
    framestyle     --> :box

    imap.data'
end

# format pixel index with its sky coordinate. ra: ind=1, dec: ind=2
function pix2sky_formatter(x, imap::Enmap; ind=1, digits=2)
    res = round(rad2deg(pix2sky(imap, float.([x, x]))[ind]), digits=digits)
    string(res)*"°"
end

end  # module
