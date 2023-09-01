@testset "Vector field                          " begin
    # generate grid
    g = MyGrid([(-1 + 0.2*i, 2π/4*j) for i in 0:3, j in 0:3])

    global u1, u2, u3

    @test_nowarn u1 = VectorField(MyScalarField(g), MyScalarField(g))
    @test_nowarn u2 = VectorField(MyScalarField, g)
    @test_nowarn u3 = VectorField(MyScalarField, g, (i, j)->i*j, (i, j)->i^2+j)

    @test grid(u1) === g

    @test parent(u2) == [u2...]
    @test u3[1] == parent(u3[1])
    u2[3] = u1[1]; @test u2[3] == u1[1]
    @test size(u2) == (3,)
    @test length(u3) == 2

    @test similar(u2) isa VectorField
    @test size(similar(u2)) == (3,)
    @test copy(u3) == u3

    # construct laplacian fields
    lap1 = MyScalarField(g)
    lap2 = MyScalarField(g)
    D = DiffMatrix([-1 + 0.2*i for i in 0:3], 3, 2)
    mul!(lap1, D, u1[1]); mul!(lap2, D, u1[2])

    # TODO: check broadcasting of functions

    # construct time derivative field
    dudt1 = MyScalarField(g)
    dudt2 = MyScalarField(g)
    dudt3 = MyScalarField(g)
    dudt1 .= one(Float64); dudt2 .= one(Float64); dudt3 .= one(Float64)

    @test laplacian!(similar(u1), u1) == VectorField(lap1, lap2)
    @test ddt!(similar(u2), u2) == VectorField(dudt1, dudt2, dudt3)
    @test cross!(similar(u2), [0, 0, 1], u2) == VectorField(-u2[2], u2[1], (v = MyScalarField(g); v .= 0.0; v))

    @test dot(u3, u1) == dot(u3[1], u1[1]) + dot(u3[2], u1[2])
    @test norm(u3) == sqrt(norm(u3[1])^2 + norm(u3[2])^2)
end
