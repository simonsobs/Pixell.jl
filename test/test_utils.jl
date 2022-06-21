@testset "dplanck" begin
    @test dplanck(98e9) == 231581854.6149218
    @test dplanck(150e9) == 398477703.66164714
end