
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
    for copy_op in (copy, deepcopy, similar)  # these should all do the same thing: NEVER keep the same WCS
        shape0, wcs0 = fullsky_geometry(π/180)
        m = Enmap(rand(shape0...), wcs0)
        m2 = copy_op(m)
        @test !(m.wcs === m2.wcs)
        m2.wcs.cdelt = [99., 99.]
        @test !(m.wcs.cdelt ≈ [99., 99.])

        m = Enmap(rand(shape0...), wcs0)
        m2 = m.^2
        @test !(m.wcs === m2.wcs)
        m2.wcs.cdelt = [99., 99.]
        @test !(m.wcs.cdelt ≈ [99., 99.])

        m = Enmap(rand(shape0...), wcs0)
        m2 = copy_op(m)
        m2 .= m
        @test !(m.wcs === m2.wcs)
        m2.wcs.cdelt = [99., 99.]  # also make sure sub-arrays aren't shared
        @test !(m.wcs.cdelt ≈ [99., 99.])
    end
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
