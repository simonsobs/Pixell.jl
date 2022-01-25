using Pixell
using Test


@testset "Enmap geometry" begin
    shape, wcs = fullsky_geometry(deg2rad(1 / 60))
    @test wcs.cdelt ≈ [-0.016666666666666666, 0.016666666666666666]
    @test wcs.crpix ≈ [10800.5, 5401.0]
    @test wcs.crval ≈ [0.008333333333333333, 0.0]

    shape, wcs = fullsky_geometry(deg2rad(1 / 61))
    @test wcs.cdelt ≈ [-0.01639344262295082, 0.01639344262295082]
    @test wcs.crpix ≈ [10980.5, 5491.0]
    @test wcs.crval ≈ [0.00819672131147541, 0.0]

    shape, wcs = fullsky_geometry(deg2rad(5); dims=(3,))
    @test shape == (72, 37, 3)
end

##
@testset "Enmap broadcasting" begin
    shape, wcs = fullsky_geometry(deg2rad(1); dims=(3,))
    A, B = rand(shape...), rand(shape...)
    ma = Enmap(A, wcs)
    mb = Enmap(B, wcs)
    @test A .+ B == ma .+ mb
    @test A .+ B == ma .+ B
    @test A .+ B == A .+ mb
    @test A .+ B .* sin.(A.^2) == (ma .+ mb .* sin.(ma.^2))
end

##
@testset "Enmap sky2pix and pix2sky" begin
    shape, wcs = fullsky_geometry(deg2rad(1))
    m = Enmap(rand(shape...), wcs)
    @test [180., -90.] ≈ pix2sky(m, [1.0, 1.0])
    @test [1.0, 0.0] ≈ sky2pix(m, pix2sky(m, [1.0, 0.0]))
    @test [13., 7.] ≈ sky2pix(m, pix2sky(m, [13., 7.]))
end

## 
@testset "nonallocating WCS info utilities" begin
    shape, wcs = fullsky_geometry(deg2rad(1))
    @test all(Pixell.crpix(wcs) .== wcs.crpix)
    @test all(Pixell.crval(wcs) .== wcs.crval)
    @test all(Pixell.cdelt(wcs) .== wcs.cdelt)
end
