@testset "dplanck" begin
    @test dplanck(98e9) == 231581854 atol = 100
    @test dplanck(150e9) == 398477703 atol = 100
end