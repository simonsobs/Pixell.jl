
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
