@testset "Abstract scalar field                 " begin
    # construct my scalar field
    g = MyGrid([(-1 + 0.2*i, 2π/4*j) for i in 0:3, j in 0:3])
    u = MyScalarField(g)

    @test MyScalarField(g, (i, j)->i*j) ≈ [0 -1π/2  -π    -3π/2;
                                            0 -2π/5  -4π/5 -6π/5;
                                            0 -3π/10 -3π/5 -9π/10;
                                            0 -π/5   -2π/5 -3π/5]
    @test size(u) == (4, 4)
    @test similar(u) isa MyScalarField && eltype(similar(u)) == Float64
    @test copy(u) == u
    @test_nowarn u[1, 2]
    @test_nowarn u[2, 4] = 2.0
    @test norm(u) == sum(u[i]^2 for i in eachindex(u))
    @test_nowarn u .= 5.234
    @test_nowarn round.(u)
end
