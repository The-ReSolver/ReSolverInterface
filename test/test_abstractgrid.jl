@testset "Abstract grid                         " begin
    # construct my grid
    grid = MyGrid([(-1 + 0.2*i, 2π/4*j) for i in 0:3, j in 0:3])

    @test points(grid) isa Matrix{Tuple{Float64, Float64}}
    @test grid == MyGrid(points(grid))
    @test !(grid == MyGrid(points(grid)[1:end-1, :]))
end
