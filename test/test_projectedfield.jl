@testset "Projected field               " begin
    modes = rand(ComplexF64, 2, 5, 2)
    a1 = @test_nowarn ProjectedField(MyGrid(), modes)
    parent(a1) .= rand(ComplexF64, 5, 2)
    @test similar(a1) isa ProjectedField{MyGrid, 2, ComplexF64, Array{ComplexF64, 3}}
    @test_nowarn norm(a1)
end
