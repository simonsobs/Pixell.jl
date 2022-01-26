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
    @test [180., -90.] ≈ collect(pix2sky(m, [1.0, 1.0]))
    @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0])))
    @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.])))

    # check that our custom implementations 
    pixcoords = 100 .* rand(2, 1024)
    skycoords = pix2sky(m, pixcoords)
    @test skycoords ≈ Pixell.WCS.pix_to_world(Pixell.getwcs(m), pixcoords)
    skycoords .= 0.0
    pix2sky!(m, pixcoords, skycoords)
    @test skycoords ≈ Pixell.WCS.pix_to_world(Pixell.getwcs(m), pixcoords)
    
    skycoords = 100 .* rand(2, 1024)
    pixcoords = sky2pix(m, skycoords)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(Pixell.getwcs(m), skycoords)
    pixcoords .= 0.0
    sky2pix!(m, skycoords, pixcoords)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(Pixell.getwcs(m), skycoords)
end

## 
@testset "nonallocating WCS info utilities" begin
    shape, wcs = fullsky_geometry(deg2rad(1))
    @test all(Pixell.crpix(wcs) .== wcs.crpix)
    @test all(Pixell.crval(wcs) .== wcs.crval)
    @test all(Pixell.cdelt(wcs) .== wcs.cdelt)
end

## 
@testset "slice_geometry" begin
    shape0, wcs0 = fullsky_geometry(deg2rad(1))
    shape, wcs = slice_geometry(shape0, wcs0, 1:4, 11:-1:2)
    @test (3, 9) == shape
    @test [-1.0, -1.0] ≈ wcs.cdelt
    @test [180.5, -79.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    
    shape, wcs = slice_geometry(shape0, wcs0, 2:shape0[1], 6:12)
    @test (358, 6) == shape
    @test [-1.0, 1.0] ≈ wcs.cdelt
    @test [179.5, 86.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 2:2:shape0[1], 1:12)
    @test (179, 11) == shape
    @test [-2.0, 1.0] ≈ wcs.cdelt
    @test [90.0, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 23:-4:3, 1:3:31)
    @test (5, 10) == shape
    @test [4.0, 3.0] ≈ wcs.cdelt
    @test [-38.75, 30.666666666666668] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
end

##
@testset "Enmap I/O" begin
    imap = read_map("data/test.fits")
    @test size(imap) == (100, 100, 3)
    @test imap.wcs.naxis == 2
    @test imap.wcs.cdelt == [-1, 1]
    @test imap.wcs.crval == [0.5, 0.0]
    @test sum(imap) ≈ 14967.2985
    # read with sel
    imap = read_map("data/test.fits", sel=(11:20,21:40,1:2))
    @test size(imap) == (10,20,2)
    # todo: add tests on IAU conversion w/ and w/o sel
end
