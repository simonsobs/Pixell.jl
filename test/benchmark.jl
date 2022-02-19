using BenchmarkTools

##

function benchmark_broadcasting()
    shape, wcs = fullsky_geometry(deg2rad(1); dims=(3,))
    A, B = rand(shape...), rand(shape...)
    ma = Enmap(A, wcs)
    mb = Enmap(B, wcs)

    @btime $mb .* $ma .* exp.($ma.^2)
    @btime $B .* $A .* exp.($A.^2)

    return (mb .* ma .* exp.(ma.^2)) == parent(mb .* ma .* exp.(ma.^2))

end

benchmark_broadcasting()

##
using Pixell, BenchmarkTools

function benchmark_sky2pix()
    shape, wcs = fullsky_geometry(deg2rad(1); dims=(3,))
    A, B = rand(shape...), rand(shape...)
    ma = Enmap(A, wcs)
    mb = Enmap(B, wcs)

    @btime sky2pix($ma, 5, 6)

end

benchmark_sky2pix()

