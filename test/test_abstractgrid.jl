@testset "Abstract grid                         " begin
    # construct my grid
    grid = MyGrid(collect(-1:0.2:1))

    @test points(grid) == [-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    @test grid == MyGrid([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    @test !(grid == MyGrid([-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]))
end
