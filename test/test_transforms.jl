
# write some simple stuff to the map
function gen_spin0!(m, powerlaw_exponent=2)
    for i in axes(m, 1)
        for j in axes(m, 2)
            m[i,j] = ( ((j-1) * size(m, 1)) + i )^powerlaw_exponent
        end
    end
end

@testset "spin-0 SHTs" begin
    shape, wcs = fullsky_geometry(10.0 * Pixell.degree)
    m = Enmap(zeros(shape), wcs)
    gen_spin0!(m)

    alms = map2alm(m; lmax=18)
    data = readdlm("data/simple_analytic_sht.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)

    @test length(map2alm(m).alm) ==5995
    
    alms = map2alm(m[6:end-2, 5:end-3]; lmax=18)
    data = readdlm("data/simple_analytic_sht_sliced.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)

    box = [10   -10;           # RA
           -5     5] * Pixell.degree  # DEC
    shape, wcs = geometry(CarClenshawCurtis, box, 1.0 * Pixell.degree)
    m = Enmap(zeros(shape), wcs)
    gen_spin0!(m, 2.5)
    alms = map2alm(m; lmax=100)
    data = readdlm("data/simple_box_analytic_sht.txt")
    ref_alms = data[:,1] + 1im .* data[:,2]
    @test (alms.alm ≈ ref_alms)
end

function gen_spin2!(m)
    for i in axes(m, 1)
        for j in axes(m, 2)
            m[i,j,1] = i^2 * j
            m[i,j,2] = i * j
        end
    end
end

@testset "spin-2 SHTs" begin
    shape, wcs = fullsky_geometry(10.0 * Pixell.degree; dims=(2,))
    m_QU = Enmap(zeros(shape), wcs)
    gen_spin2!(m_QU)
    alm_e, alm_b = map2alm(m_QU; lmax=3 * 36)
    data = readdlm("data/simple_pol_analytic_sht.txt")
    ref_alm_e = data[:,1] + 1im .* data[:,2]
    ref_alm_b = data[:,3] + 1im .* data[:,4]
    @test (ref_alm_e ≈ alm_e.alm)
    @test (ref_alm_b ≈ alm_b.alm)

    alm_e, alm_b = map2alm((m_QU[:,:,1], m_QU[:,:,2]); lmax=3 * 36)
    @test (ref_alm_e ≈ alm_e.alm)
    @test (ref_alm_b ≈ alm_b.alm)

    shape, wcs = fullsky_geometry(10.0 * Pixell.degree; dims=(3,))
    m_IQU = Enmap(zeros(shape), wcs)
    gen_spin0!(@view(m_IQU[:,:,1]))
    gen_spin2!(@view(m_IQU[:,:,2:3]))
    alm_t, alm_e, alm_b = map2alm(m_IQU; lmax=3 * 36)
    
    data = readdlm("data/simple_analytic_sht_fullalm.txt")
    ref_alm_t = data[:,1] + 1im .* data[:,2]
    @test (ref_alm_t ≈ alm_t.alm)
    @test (ref_alm_e ≈ alm_e.alm)
    @test (ref_alm_b ≈ alm_b.alm)
    
    alm_t, alm_e, alm_b = map2alm((m_IQU[:,:,1], m_IQU[:,:,2], m_IQU[:,:,3]); lmax=3 * 36)
end
