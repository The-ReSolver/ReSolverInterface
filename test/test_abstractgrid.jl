@testset "Grid                          " begin
    g = MyGrid()
    pts = @test_nowarn points(g)
    @test g == MyGrid()
end
