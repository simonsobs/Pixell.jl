using Pixell
using Test
import Pixell: degree, arcminute

@testset "Enmap geometry" begin
    for W in (Pixell.WCS.WCSTransform, CarClenshawCurtis{Float64})
        shape, wcs = fullsky_geometry(W, deg2rad(1 / 60))
        @test collect(wcs.cdelt) ≈ [-0.016666666666666666, 0.016666666666666666]
        @test collect(wcs.crpix) ≈ [10800.5, 5401.0]
        @test collect(wcs.crval) ≈ [0.008333333333333333, 0.0]

        shape, wcs = fullsky_geometry(W, deg2rad(1 / 61))
        @test collect(wcs.cdelt) ≈ [-0.01639344262295082, 0.01639344262295082]
        @test collect(wcs.crpix) ≈ [10980.5, 5491.0]
        @test collect(wcs.crval) ≈ [0.00819672131147541, 0.0]

        shape, wcs = fullsky_geometry(W, deg2rad(5); dims=(3,))
        @test shape == (72, 37, 3)

        box = [10  -10;           # RA
               -5    5] * degree  # DEC
        shape, wcs = geometry(W, box, 0.5 * arcminute)
        @test shape == (2400, 1200)
        @test collect(wcs.cdelt) ≈ [-0.008333333333333333,  0.008333333333333333]
        @test collect(wcs.crpix) == [1201, 601]
        @test collect(wcs.crval) == [0.0, 0.0]

        box = [11  -10;           # RA
               -6    5] * degree  # DEC
        shape, wcs = geometry(W, box, 0.5 * arcminute)
        @test shape == (2520, 1320)
        @test collect(wcs.cdelt) ≈ [-0.008333333333333333,  0.008333333333333333]
        @test collect(wcs.crpix) ≈ [1261, 721]
        @test collect(wcs.crval) ≈ [0.5, 0.0]

        box = [10   -4;           # RA
               -3    5] * degree  # DEC
        shape, wcs = geometry(W, box, 0.5 * arcminute)
        @test shape == (1680, 960)
        @test collect(wcs.crpix) ≈ [841, 361]
        @test collect(wcs.crval) ≈ [3.0, 0.0]
    end
end

##
wrap(ra_dec_vec) = [mod(ra_dec_vec[1], 2π), mod(ra_dec_vec[2], π)]
@testset "Enmap sky2pix and pix2sky" begin
    for W in (Pixell.WCS.WCSTransform, CarClenshawCurtis{Float64})
        shape, wcs = fullsky_geometry(W, deg2rad(1))
        m = Enmap(rand(shape...), wcs)
        # in this test, wrap to angles in [0, 2π] and [0, π] for RA and DEC
        @test [3.12413936, -1.55334303] ≈ collect(pix2sky(m, [2.0, 2.0]))
        @test [2.96705973, -1.79768913] ≈ collect(pix2sky(m, [11.0, -12.0]))
        @test [2.44346095, -2.0943951] ≈ collect(pix2sky(m, [41.0, -29.0]))
        @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0])))
        @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.])))
        
        # test passing a vector
        @test [[1.0], [0.0]] ≈ collect(sky2pix(m, pix2sky(m, [1.0], [0.0])...))
        @test [[13.], [7.]] ≈ collect(sky2pix(m, pix2sky(m, [13.], [7.])...))
        
        # try some crazy angles
        @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0]) .+ (2π, 6π)))
        @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.]) .+ (12π, 16π)))

        # check our custom implementations
        pixcoords = π .* rand(2, 1024)
        skycoords = pix2sky(m, pixcoords; safe=false)
        shape, wcs_generic = fullsky_geometry(Pixell.WCS.WCSTransform, deg2rad(1))
        @test skycoords ≈ Pixell.WCS.pix_to_world(wcs_generic, pixcoords) .* (π/180)
        skycoords .= 0.0
        pix2sky!(m, pixcoords, skycoords; safe=false)
        @test skycoords ≈ Pixell.WCS.pix_to_world(wcs_generic, pixcoords) .* (π/180)
        
        skycoords = π .* rand(2, 1024)
        pixcoords = sky2pix(m, skycoords; safe=false)
        @test pixcoords ≈ Pixell.WCS.world_to_pix(wcs_generic, skycoords .* (180/π))
        pixcoords .= 0.0
        sky2pix!(m, skycoords, pixcoords; safe=false)
        @test pixcoords ≈ Pixell.WCS.world_to_pix(wcs_generic, skycoords .* (180/π))

        box = [10   -10;           # RA
            -5     5] * degree  # DEC
        shape, wcs = geometry(W, box, 1 * degree)
        m = Enmap(ones(shape), wcs)
        @test [sky2pix(m, deg2rad(ra), 0.0; safe=true)[1]
            for ra in 178:182] ≈ [-167., -168., -169.,  190.,  189.]
    end
end

##
@testset "Enmap sky2pix and pix2sky WCSlib fallbacks" begin
    shape, wcs = fullsky_geometry(Pixell.WCS.WCSTransform, deg2rad(1))
    m = Enmap(rand(shape...), wcs)
    # in this test, wrap to angles in [0, 2π] and [0, π] for RA and DEC
    @test [3.12413936, -1.55334303] ≈ collect(pix2sky(m, [2.0, 2.0]))
    @test [2.96705973, -1.79768913] ≈ collect(pix2sky(m, [11.0, -12.0]))
    @test [2.44346095, -2.0943951] ≈ collect(pix2sky(m, [41.0, -29.0]))
    @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0])))
    @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.])))

    
    @test [1.0, 0.0] ≈ collect(sky2pix(m, pix2sky(m, [1.0, 0.0]) .+ (2π, 6π)))
    @test [13., 7.] ≈ collect(sky2pix(m, pix2sky(m, [13., 7.]) .+ (12π, 16π)))

    # check our custom implementations
    pixcoords = π .* rand(2, 1024)
    skycoords = pix2sky(m, pixcoords; safe=false)
    shape, wcs_generic = fullsky_geometry(Pixell.WCS.WCSTransform, deg2rad(1))
    @test skycoords ≈ Pixell.WCS.pix_to_world(wcs_generic, pixcoords) .* (π/180)
    skycoords .= 0.0
    pix2sky!(m, pixcoords, skycoords; safe=false)
    @test skycoords ≈ Pixell.WCS.pix_to_world(wcs_generic, pixcoords) .* (π/180)
    
    skycoords = π .* rand(2, 1024)
    pixcoords = sky2pix(m, skycoords; safe=false)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(wcs_generic, skycoords .* (180/π))
    pixcoords .= 0.0
    sky2pix!(m, skycoords, pixcoords; safe=false)
    @test pixcoords ≈ Pixell.WCS.world_to_pix(wcs_generic, skycoords .* (180/π))

    box = [10   -10;           # RA
           -5     5] * degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis, box, 1 * degree)
    m = Enmap(ones(shape), wcs)
    @test [sky2pix(m, deg2rad(ra), 0.0; safe=true)[1]
        for ra in 178:182] ≈ [-167., -168., -169.,  190.,  189.]

end

##
@testset "sky2pix every single pixel" begin
    box = [10   -10;           # RA
        -5     5] * degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis, box, 1 * degree)
    m = Enmap(ones(shape), wcs)
    for i in 1:shape[1]
        for j in 1:shape[2]
            ra, dec = pix2sky(m, i, j)
            ra_unsafe, dec_unsafe = pix2sky(m, i, j; safe=false)
            @test ra ≈ ra_unsafe 
            @test dec ≈ dec_unsafe
            
            ra_pix, dec_pix = sky2pix(m, ra, dec)
            ra_pix_unsafe, dec_pix_unsafe = sky2pix(m, ra, dec; safe=false)

            @test ra_pix ≈ ra_pix_unsafe
            @test dec_pix ≈ dec_pix_unsafe
        end
    end

    box = [179   -179;           # RA
           -89     89] * degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis, box, 1 * degree)
    m = Enmap(ones(shape), wcs)
    for i in 1:shape[1]
        for j in 1:shape[2]
            ra, dec = pix2sky(m, i, j)
            ra_unsafe, dec_unsafe = pix2sky(m, i, j; safe=false)
            @test ra ≈ ra_unsafe 
            @test dec ≈ dec_unsafe
            @test -π ≤ ra ≤ π
            @test -π/2 ≤ dec ≤ π/2
            
            ra_pix, dec_pix = sky2pix(m, ra, dec)
            ra_pix_unsafe, dec_pix_unsafe = sky2pix(m, ra, dec; safe=false)

            @test ra_pix ≈ ra_pix_unsafe
            @test dec_pix ≈ dec_pix_unsafe
            @test 1 ≤ ra_pix ≤ shape[1]
            @test 1 ≤ dec_pix ≤ shape[2]
        end
    end
    
    # determine if safe angles are on the sky
    shape, wcs = fullsky_geometry(deg2rad(1))
    m = Enmap(ones(shape), wcs)
    for i in 1:shape[1]
        for j in 1:shape[2]
            ra, dec = pix2sky(m, i, j)
            ra_unsafe, dec_unsafe = pix2sky(m, i, j; safe=false)
            @test -π ≤ ra ≤ π
            @test -π/2 ≤ dec ≤ π/2

            ra_pix, dec_pix = sky2pix(m, ra, dec)
            @test 1 ≤ ra_pix ≤ shape[1]
            @test 1 ≤ dec_pix ≤ shape[2]
            ra_pix, dec_pix = sky2pix(m, ra_unsafe, dec_unsafe)
            @test 1 ≤ ra_pix ≤ shape[1]
            @test 1 ≤ dec_pix ≤ shape[2]
        end
    end
    
end

## 
@testset "nonallocating WCS info utilities" begin
    shape, wcs = fullsky_geometry(deg2rad(1))
    @test all(Pixell.getcrpix(wcs) .== collect(wcs.crpix))
    @test all(Pixell.getcrval(wcs) .== collect(wcs.crval))
    @test all(Pixell.getcdelt(wcs) .== collect(wcs.cdelt))
end

## 
@testset "slice_geometry" begin
    for W in (Pixell.WCS.WCSTransform, CarClenshawCurtis{Float64})
        shape0, wcs0 = fullsky_geometry(W, deg2rad(1))
        shape, wcs = slice_geometry(shape0, wcs0, 1:3, 11:-1:3)
        @test (3, 9) == shape
        @test [-1.0, -1.0] ≈ collect(wcs.cdelt)
        @test [180.5, -79.0] ≈ collect(wcs.crpix)
        @test [0.5, 0.0] ≈ collect(wcs.crval)

        shape, wcs = slice_geometry(shape0, wcs0, 2:shape0[1], 6:11)
        @test (359, 6) == shape
        @test [-1.0, 1.0] ≈ collect(wcs.cdelt)
        @test [179.5, 86.0] ≈ collect(wcs.crpix)
        @test [0.5, 0.0] ≈ collect(wcs.crval)

        shape, wcs = slice_geometry(shape0, wcs0, 2:2:(shape0[1]-1),1:11)
        @test (179, 11) == shape
        @test [-2.0, 1.0] ≈ collect(wcs.cdelt)
        @test [90.0, 91.0] ≈ collect(wcs.crpix)
        @test [0.5, 0.0] ≈ collect(wcs.crval)

        shape, wcs = slice_geometry(shape0, wcs0, 23:-4:6, 1:3:28)
        @test (5, 10) == shape
        @test [4.0, 3.0] ≈ collect(wcs.cdelt)
        @test [-38.75, 30.666666666666668] ≈ collect(wcs.crpix)
        @test [0.5, 0.0] ≈ collect(wcs.crval)

        shape, wcs = slice_geometry(shape0, wcs0, 3:3, 1:3:28)
        @test (1, 10) == shape
        @test [-1.0, 3.0] ≈ collect(wcs.cdelt)
        @test [178.5, 30.666666666666668] ≈ collect(wcs.crpix)
        @test [0.5, 0.0] ≈ collect(wcs.crval)
    end
end

## 
@testset "pixel skyareas" begin
    for W in (Pixell.WCS.WCSTransform, CarClenshawCurtis{Float64})
        shape, wcs = fullsky_geometry(W, deg2rad(1))
        @test skyarea(shape, wcs) ≈ 4π

        shape′, wcs′ = slice_geometry(shape, wcs, 9:-1:3, 1:2)
        @test skyarea(shape′, wcs′) ≈ 4.1865652086145036e-05

        box = [10   -10;           # RA
               -5     5] * degree  # DEC
        shape, wcs = geometry(CarClenshawCurtis, box, (1/60) * degree)
        @test skyarea(shape, wcs) ≈ 0.06084618627514243
    end
end

@testset "pixareamap" begin
    shape, wcs = fullsky_geometry(π/180)
    pm = pixareamap(shape, wcs)
    pm_py = readdlm("data/fullsky_pixareas.dat")
    @test sum(abs.(pm[1,:] .- pm_py)) < 100eps()

    box = [10  -10;           # RA
           -5    5] * Pixell.degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 5. * Pixell.arcminute)
    pm = pixareamap(shape, wcs)
    pm_py = readdlm("data/box_pixareas.dat")
    @test sum(abs.(pm[1,:] .- pm_py)) < 100eps()
end
