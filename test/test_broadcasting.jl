@testset "Scalar field broacasting      " begin
    u = MyField(MyGrid(), rand(ComplexF64, 2, 2))
    v = MyField(MyGrid(), parent(u) .+ Complex(1.0, 2.0))
    @test u .+ Complex(1.0, 2.0) == v
    @test copy(u) == u
end

@testset "Vector field broadcasting     " begin
    g = MyGrid()
    u1 = VectorField(MyField(g, ones(2, 2)), MyField(g, ones(2, 2)))
    u2 = VectorField(MyField, g)
    u3 = VectorField(MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)), MyField(g, rand(ComplexF64, 2, 2)))

    @test u1 .+ (1.0+2.0im) == VectorField(MyField(g, [2+2im 2+2im; 2+2im 2+2im]), MyField(g, [2+2im 2+2im; 2+2im 2+2im]))
    @test u1 .+ [2+1im, 5] == VectorField(MyField(g, [3+1im 3+1im; 3+1im 3+1im]), MyField(g, [6+0im 6+0im; 6+0im 6+0im]))
    @test u1 .+ VectorField(MyField(g, 2*ones(2, 2)), MyField(g, 3im*ones(2, 2))) == VectorField(MyField(g, [3+0im 3+0im; 3+0im 3+0im]), MyField(g, [1+3im 1+3im; 1+3im 1+3im]))
end
