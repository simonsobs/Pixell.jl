
@testset "different distance transform implementations" begin

    box = [20   -20;           # RA
          -10     10] * degree  # DEC
    shape, wcs = geometry(Pixell.WCS.WCSTransform, box, (1/2) * degree)

    for kk in 1:300
        m = Enmap(ones(shape), wcs)
        # m[2,2] = 0.0
        for i in 1:30
            m[rand(2:(size(m,1)-1)), rand(2:(size(m,2)-1))] = 0.0
        end
        
        distmap = distance_transform(ExactSeqSDT(), m)
        distmap_bf = distance_transform(BruteForceSDT(), m)
        distmap_approx = distance_transform(ApproxSeqSDT(), m)

        @test sum(abs.((distmap_bf .- distmap))) < 1e-13
        @test (sum(distmap_bf .!= distmap_approx) / prod(size(distmap))) < 0.2
    end

end

##
@testset "distance transform metric" begin
    box = [20   -20;           # RA
           0     10] * degree  # DEC
    shape, wcs = geometry(Pixell.WCS.WCSTransform, box, (1/2) * degree)

    m = Enmap(ones(shape), wcs)
    m[1,1] = 0.0
    distmap = distance_transform(ExactSeqSDT(), m)
    
    αs = pix2sky(m, collect(1:size(m,1)), ones(size(m,1)))[1]
    δs = pix2sky(m, ones(size(m,2)), collect(1:size(m,2)))[2]

    for i in axes(m,1)
        @test (αs[1] - αs[i]) ≈ distmap[i,1]
    end
    for j in axes(m,2)
        @test (δs[j] - δs[1]) ≈ distmap[1,j]
    end
end
