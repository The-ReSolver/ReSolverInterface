@testset "Projected field               " begin
    g = MyGrid()
    modes = rand(ComplexF64, 5, 5, 5)
    a1 = @test_nowarn ProjectedField(MyField(g), modes)
    @test similar(a1) isa ProjectedField{MyField, 2, Float64, Array{ComplexF64, 3}}
end
