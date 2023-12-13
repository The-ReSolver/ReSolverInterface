@testset "Scalar field                  " begin
    # u = MyField(MyGrid(), rand(ComplexF64, 2, 2))
    u = MyField(MyGrid(), ones(ComplexF64, 2, 2))
    v = MyField(MyGrid(), parent(u) .+ Complex(1.0, 2.0))
    @test size(u) == (2, 2)
    @test IndexStyle(MyField) == Base.IndexLinear()
    @test u .+ Complex(1.0, 2.0) == v
    @test copy(u) == u
    @test u[1] == u.data[1]
    @test u[2, 1] == u.data[2, 1]
    @test_nowarn u[2, 2] = 2.0
    @test u[2, 2] == 2.0
    @test norm(u) == sqrt(sum(dot(ui, ui) for ui in u))
end
