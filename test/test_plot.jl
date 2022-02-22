@testset "Enmap plot" begin
    using Plots; gr()
    using Pixell.Enplot
    register_colors!()
    imap = read_map("data/test.fits", sel=(:,:,1))
    @test (plot(imap, color=:planck); true)
end