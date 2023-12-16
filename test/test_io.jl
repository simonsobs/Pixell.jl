

@testset "Enmap I/O" begin
    for trim in (true, false)
        imap = read_map("data/test.fits"; trim=trim)
        @test size(imap) == (100, 100, 3)
        @test imap.wcs.naxis == 2
        @test collect(imap.wcs.cdelt) == [-1, 1]
        @test collect(imap.wcs.crval) == [0.5, 0.0]
        @test sum(imap) ≈ 14967.2985
        # read with sel
        imap = read_map("data/test.fits", sel=(11:20,21:40,1:2); trim=trim)
        @test size(imap) == (10,20,2)
    end
    # todo: add tests on IAU conversion w/ and w/o sel
end

##
@testset "Enmap ClenshawCurtis vs Fejer1" begin
    
    imap1 = read_map("data/cc.fits"; trim=false)
    @test Pixell.isfejer1(imap1.wcs) == false
    @test Pixell.isclenshawcurtis(imap1.wcs)

    imap2= read_map("data/fejer1.fits"; trim=false)
    @test Pixell.isfejer1(imap2.wcs)
    @test Pixell.isclenshawcurtis(imap2.wcs) == false

    imap1 = read_map("data/cc.fits"; trim=true)
    @test Pixell.isfejer1(imap1.wcs) == false
    @test Pixell.isclenshawcurtis(imap1.wcs)

    imap2= read_map("data/fejer1.fits"; trim=true)
    @test Pixell.isfejer1(imap2.wcs)
    @test Pixell.isclenshawcurtis(imap2.wcs) == false

end
