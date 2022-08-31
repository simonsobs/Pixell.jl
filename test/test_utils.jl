# using Test, DelimitedFiles, LinearAlgebra, Pixell

@testset "dplanck" begin
    @test dplanck(98e9) ≈ 231581854 atol = 100
    @test dplanck(150e9) ≈ 398477703 atol = 100
end

@testset "FFTLog" begin
    N = 64
    μ = 0
    q = 0.0
    r₀ = 1.0
    L = 8.0
    Nhalf = N ÷ 2
    n = range(-Nhalf,Nhalf,length=N)
    r = r₀ .* 10 .^ (n .* L ./ N )
    pl = Pixell.plan_fftlog(r, μ, q, 1.0; kropt=true)
    aₙ = r .^ (μ + 1) .* exp.(-r.^2 / 2)
    y = similar(r, ComplexF64)
    fftdata = readdlm("data/fftlog_example.txt", ' ', Float64, '\n')

    # test forward
    mul!(y, pl, aₙ)
    f_ref = fftdata[:,2]
    @test all(abs.(y .- f_ref) .< 1e-15)
    @test isapprox(y, f_ref)

    # test backward
    y2 = similar(r, ComplexF64)
    ldiv!(y2, pl, y)
    @test all(abs.(y2 .- aₙ) .< 1e-15)
end

@testset begin
    rft = RadialFourierTransform(n=256, pad=128)
    rftdata = readdlm("data/radialfouriertransform.txt")

    h = real2harm(rft, r -> 1/r)
    @test maximum(abs.(1 .- h ./ rftdata[:,1])) < 1000eps()

    h = harm2real(rft, r -> 1/r)
    @test maximum(abs.(1 .- h ./ rftdata[:,2])) < 1000eps()

    
    h = real2harm(rft, 1 ./ rft.r)
    @test maximum(abs.(1 .- h ./ rftdata[:,1])) < 1000eps()

    h = harm2real(rft, 1 ./ rft.revl)
    @test maximum(abs.(1 .- h ./ rftdata[:,2])) < 1000eps()

end
