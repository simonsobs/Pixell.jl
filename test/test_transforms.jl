
# write some simple stuff to the map
function gen_simple!(m, powerlaw_exponent)
    for i in axes(m, 1)
        for j in axes(m, 2)
            m[i,j] = ( ((j-1) * size(m, 1)) + i )^powerlaw_exponent
        end
    end
end

@testset "Scalar SHTs" begin
    shape, wcs = fullsky_geometry(10.0 * Pixell.degree)
    m = Enmap(zeros(shape), wcs)
    gen_simple!(m, 2)

    alms = map2alm(m; lmax=18)
    data = readdlm("data/simple_analytic_sht.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)
    
    alms = map2alm(m[6:end-2, 5:end-3]; lmax=18)
    data = readdlm("data/simple_analytic_sht_sliced.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)

    
    box = [10   -10;           # RA
           -5     5] * degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis, box, 1.0 * degree)
    m = Enmap(zeros(shape), wcs)
    gen_simple!(m, 2.5)
    alms = map2alm(m; lmax=100)
    data = readdlm("data/simple_box_analytic_sht.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)
end
