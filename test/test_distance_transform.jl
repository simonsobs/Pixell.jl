
@testset "fast exact vs bruteforce distance transform" begin

        
    box = [20   -20;           # RA
          -10     10] * degree  # DEC
    shape, wcs = geometry(Pixell.WCS.WCSTransform, box, (1/2) * degree)

    for kk in 1:300
        m = Enmap(ones(shape), wcs)
        # m[2,2] = 0.0
        for i in 1:30
            m[rand(2:(size(m,1)-1)), rand(2:(size(m,2)-1))] = 0.0
        end
        
        distmap = distance_transform(ExactSeqSDT(1), m)
        distmap_bf = distance_transform(BruteForceSDT(), m)
        distmap_approx = distance_transform(ApproxSeqSDT(), m)

        @test sum(abs.((distmap_bf .- distmap))) < 1e-13
        @test (sum(distmap_bf .!= distmap_approx) / prod(size(distmap))) < 0.1
    end

end

