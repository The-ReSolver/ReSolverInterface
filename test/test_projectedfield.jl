@testset "Projected field               " begin
    modes = rand(ComplexF64, 5, 5, 5)
    a1 = @test_nowarn ProjectedField(modes, rand(ComplexF64, 5, 2))
    @test similar(a1) isa ProjectedField{2, ComplexF64, Array{ComplexF64, 3}}
end
