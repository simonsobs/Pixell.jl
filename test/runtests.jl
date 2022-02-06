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

    ma .= 1.0
    @test all(ma .≈ 1.0)
    ma .= mb
    @test all(ma .≈ mb)
    ma[1,:,3] .= 2.0
    @test all(ma[1,:,3] .≈ 2.0)
    ma[:,end,3] .= 3.0
    @test all(ma[:,end,3] .≈ 3.0)
    ma[:,1,:] .= mb[:,2,:]
    @test all(ma.data[:,1,:] .≈ mb.data[:,2,:])

    A, B = rand(shape...), rand(shape...)
    ma = Enmap(A, wcs)
    mb = Enmap(B, wcs)
    mb[:,:,1] .= ma[:,:,1]
    @test all(ma.data[:,:,1] .≈ mb.data[:,:,1])
    @test !all(ma.data[:,:,2] .≈ mb.data[:,:,2])
    @test !all(ma.data[:,:,3] .≈ mb.data[:,:,3])
    mb[:,:,:] .= ma[:,:,:]
    @test all(ma.data .≈ mb.data)

    mv = @view ma[1,:,1]
    @test Pixell.getwcs(mv) == Pixell.NoWCS()
    mv = ma[1,:,1]
    @test Pixell.getwcs(mv) == Pixell.NoWCS()
    mv = ma[1,:,:]
    @test Pixell.getwcs(mv) == Pixell.NoWCS()
    mv = ma[1:5,:,1]
    @test Pixell.getwcs(mv) != Pixell.NoWCS()
    mv = ma[1:5,:,:]
    @test Pixell.getwcs(mv) != Pixell.NoWCS()
end

##

wrap(ra_dec_vec) = [mod(ra_dec_vec[1], 2π), mod(ra_dec_vec[2], π)]
@testset "Enmap sky2pix and pix2sky" begin
    shape, wcs = fullsky_geometry(deg2rad(1))
    m = Enmap(rand(shape...), wcs)
    # in this test, wrap to angles in [0, 2π] and [0, π] for RA and DEC
    @test [3.12413936, -1.55334303] ≈ collect(pix2sky(m, [2.0, 2.0]))
    @test [2.96705973, -1.79768913] ≈ collect(pix2sky(m, [11.0, -12.0]))
    @test [2.44346095, -2.0943951] ≈ collect(pix2sky(m, [41.0, -29.0]))
    @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0])))
    @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.])))

    # check that our custom implementations 
    pixcoords = π .* rand(2, 1024)
    skycoords = pix2sky(m, pixcoords; safe=false)
    @test skycoords ≈ Pixell.WCS.pix_to_world(Pixell.getwcs(m), pixcoords) .* (π/180)
    skycoords .= 0.0
    pix2sky!(m, pixcoords, skycoords; safe=false)
    @test skycoords ≈ Pixell.WCS.pix_to_world(Pixell.getwcs(m), pixcoords) .* (π/180)
    
    skycoords = π .* rand(2, 1024)
    pixcoords = sky2pix(m, skycoords; safe=false)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(Pixell.getwcs(m), skycoords .* (180/π))
    pixcoords .= 0.0
    sky2pix!(m, skycoords, pixcoords; safe=false)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(Pixell.getwcs(m), skycoords .* (180/π))
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
    shape, wcs = slice_geometry(shape0, wcs0, 1:3, 11:-1:3)
    @test (3, 9) == shape
    @test [-1.0, -1.0] ≈ wcs.cdelt
    @test [180.5, -79.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 2:shape0[1], 6:11)
    @test (359, 6) == shape
    @test [-1.0, 1.0] ≈ wcs.cdelt
    @test [179.5, 86.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 2:2:(shape0[1]-1),1:11)
    @test (179, 11) == shape
    @test [-2.0, 1.0] ≈ wcs.cdelt
    @test [90.0, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 23:-4:6, 1:3:28)
    @test (5, 10) == shape
    @test [4.0, 3.0] ≈ wcs.cdelt
    @test [-38.75, 30.666666666666668] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval

    shape, wcs = slice_geometry(shape0, wcs0, 3:3, 1:3:28)
    @test (1, 10) == shape
    @test [-1.0, 3.0] ≈ wcs.cdelt
    @test [178.5, 30.666666666666668] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
end

## 
@testset "Enmap slicing" begin
    shape0, wcs0 = fullsky_geometry(π/180)
    m = Enmap(rand(shape0...), wcs0)

    # regular slicing
    m_sliced = m[5:10, 1:end]
    wcs = getwcs(m_sliced)
    @test [-1.0, 1.0] ≈ wcs.cdelt
    @test [176.5, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[5:10, 1:end] ≈ m_sliced

    # backwards slicing
    m_sliced = m[1:12, end:-1:begin]
    wcs = getwcs(m_sliced)
    @test [-1.0, -1.0] ≈ wcs.cdelt
    @test [180.5, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[1:12, end:-1:begin] ≈ m_sliced

    # non-unit steps
    m_sliced = m[3:5:24, 39:-3:2]
    wcs = getwcs(m_sliced)
    @test [-5.0, -3.0] ≈ wcs.cdelt
    @test [36.1, -16.666666666666668] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[3:5:24, 39:-3:2] ≈ m_sliced
end

## 
@testset "Enmap view slicing" begin
    shape0, wcs0 = fullsky_geometry(π/180)
    m = Enmap(rand(shape0...), wcs0)

    # regular slicing
    m_sliced = @view m[5:10, 1:end]
    wcs = getwcs(m_sliced)
    @test [-1.0, 1.0] ≈ wcs.cdelt
    @test [176.5, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[5:10, 1:end] ≈ m_sliced

    # backwards slicing
    m_sliced = @view m[1:12, end:-1:begin]
    wcs = getwcs(m_sliced)
    @test [-1.0, -1.0] ≈ wcs.cdelt
    @test [180.5, 91.0] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[1:12, end:-1:begin] ≈ m_sliced

    # non-unit steps
    m_sliced = @view m[3:5:24, 39:-3:2]
    wcs = getwcs(m_sliced)
    @test [-5.0, -3.0] ≈ wcs.cdelt
    @test [36.1, -16.666666666666668] ≈ wcs.crpix
    @test [0.5, 0.0] ≈ wcs.crval
    @test (m.data)[3:5:24, 39:-3:2] ≈ m_sliced
end

## WCS should be not be shared under deepcopy, broadcasting, or broadcasted assignment
@testset "Enmap copying behavior" begin  
    shape0, wcs0 = fullsky_geometry(π/180)
    m = Enmap(rand(shape0...), wcs0)
    m2 = deepcopy(m)
    @test !(m.wcs === m2.wcs)
    m2.wcs.cdelt = [99., 99.]
    @test !(m.wcs.cdelt ≈ [99., 99.])

    m = Enmap(rand(shape0...), wcs0)
    m2 = m.^2
    @test !(m.wcs === m2.wcs)
    m2.wcs.cdelt = [99., 99.]
    @test !(m.wcs.cdelt ≈ [99., 99.])

    m = Enmap(rand(shape0...), wcs0)
    m2 = deepcopy(m)
    m2 .= m
    @test !(m.wcs === m2.wcs)
    m2.wcs.cdelt = [99., 99.]  # also make sure sub-arrays aren't shared
    @test !(m.wcs.cdelt ≈ [99., 99.])
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
