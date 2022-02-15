

@testset "Enmap I/O" begin
    imap = read_map("data/test.fits")
    @test size(imap) == (100, 100, 3)
    @test imap.wcs.naxis == 2
    @test imap.wcs.cdelt == [-1, 1]
    @test imap.wcs.crval == [0.5, 0.0]
    @test sum(imap) ≈ 14967.2985
    # read with sel
    imap = read_map("data/test.fits", sel=(11:20,21:40,1:2))
    @test size(imap) == (10,20,2)
    # todo: add tests on IAU conversion w/ and w/o sel
end
