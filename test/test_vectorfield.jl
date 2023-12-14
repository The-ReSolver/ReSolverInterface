@testset "Vector field                  " begin
    g = MyGrid()
    u1 = @test_nowarn VectorField(MyField(g), MyField(g))
    u2 = @test_nowarn VectorField(MyField, g)
    u3 = @test_nowarn VectorField(MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)))
    u4 = @test_nowarn VectorField(MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)))

    @test grid(u1) === g

    @test parent(u2) == [u2...]
    @test u2[1] == parent(u2[1])
    u2[3] = u1[1]; @test u2[3] == u1[1]
    @test size(u2) == (3,)
    @test length(u1) == 2

    @test similar(u2) isa VectorField{3, MyField}
    @test copy(u2) == u2

    @test cross!(similar(u2), [0, 0, 1], u2) == VectorField(-u2[2], u2[1], MyField(g))

    @test dot(u3, u4) ≈ dot(u3[1], u4[1]) + dot(u3[2], u4[2]) + dot(u3[3], u4[3])
    @test norm(u3) ≈ sqrt(norm(u3[1])^2 + norm(u3[2])^2 + norm(u3[3])^2)
end
