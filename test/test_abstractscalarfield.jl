@testset "Scalar field construction     " begin
    # construct grid
    Ny = rand(3:20)
    Nz = rand(3:20)
    Nt = rand(3:20)
    y = rand(Ny)
    D1 = rand(Ny, Ny)
    D2 = rand(Ny, Ny)
    w = rand(Ny)
    ω = abs(randn())
    β = abs(rand())
    g = RPCFGrid(y, Nz, Nt, D1, D2, w, ω, β)

    # construct field with random elements
    u = RPCFScalarField(g, rand(ComplexF64, Ny, (Nz >> 1) + 1, Nt))

    # construct field and test
    @test u isa RPCFScalarField{Ny, Nz, Nt, Float64, Matrix{Float64}}
    @test size(u) == (Ny, (Nz >> 1) + 1, Nt)
    @test similar(u) == RPCFScalarField(g)
    @test copy(u) == u
    @test zero(u) == RPCFScalarField(g, zeros(ComplexF64, Ny, (Nz >> 1) + 1, Nt))
    @test_nowarn u[rand(1:Ny), rand(1:((Nz >> 1) + 1)), rand(1:Nt)]
    @test_nowarn u[rand(1:Ny), rand(1:((Nz >> 1) + 1)), rand(1:Nt)] = 2.0
end

@testset "Scalar field interface        " begin

end

@testset "Scalar field broadcasting     " begin

end

@testset "Scalar field transform        " begin

end

@testset "Scalar field multiplication   " begin

end

@testset "Scalar field vector calculus  " begin

end

@testset "Scalar field inner product    " begin

end
