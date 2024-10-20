nalloc(a, b, c) = @allocated a .= 3.0.*b .- c./2

@testset "Scalar field broacasting      " begin
    u = MyField(MyGrid(), rand(ComplexF64, 2, 2))
    v = MyField(MyGrid(), parent(u) .+ Complex(1.0, 2.0))
    @test u .+ Complex(1.0, 2.0) == v
    @test copy(u) == u
end

@testset "Vector field broadcasting     " begin
    g = MyGrid()
    u1 = VectorField(MyField(g, ones(2, 2)), MyField(g, ones(2, 2)))
    u2 = VectorField(MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)))

    @test_nowarn u2 .= 0.0

    @test nalloc(u2, similar(u2), similar(u2)) == 0

    @test u1 .+ (1.0+2.0im) == VectorField(MyField(g, [2+2im 2+2im; 2+2im 2+2im]), MyField(g, [2+2im 2+2im; 2+2im 2+2im]))
    @test u1.*[2+1im, 5] == VectorField(MyField(g, ones(2, 2)).*[2+1im, 5], MyField(g, ones(2, 2)).*[2+1im, 5])
    @test u1 .+ VectorField(MyField(g, 2*ones(2, 2)), MyField(g, 3im*ones(2, 2))) == VectorField(MyField(g, 3.0.*ones(2, 2)), MyField(g, (1+3im).*ones(2, 2)))
    @test u1.*MyField(g, (1+2im).*ones(2, 2)) == VectorField(MyField(g, (1+2im).*ones(2, 2)), MyField(g, (1+2im).*ones(2, 2)))
end

@testset "Projected field broadcasting  " begin
    g = MyGrid()
    modes = rand(ComplexF64, 5, 5, 5)
    a1 = ProjectedField(g, modes)
    a1 .= rand(ComplexF64, 5, 2)
    a2 = ProjectedField(MyField(g), modes)
    a3 = ProjectedField(VectorField(g), modes)

    @test nalloc(a1, a2, a3) == 0

    @test a1 .+ a2 == ProjectedField(modes, MyGrid(), a1.field + a2.field)
end
