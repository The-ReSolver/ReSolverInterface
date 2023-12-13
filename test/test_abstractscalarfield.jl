@testset "Scalar field construction     " begin
    # construct field
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
    u = RPCFScalarField(g, rand(ComplexF64, Ny, (Nz >> 1) + 1, Nt))

    @test u isa RPCFScalarField{Ny, Nz, Nt, Float64, Matrix{Float64}}
    @test size(u) == (Ny, (Nz >> 1) + 1, Nt)
    @test similar(u) == RPCFScalarField(g)
    @test copy(u) == u
    @test zero(u) == RPCFScalarField(g, zeros(ComplexF64, Ny, (Nz >> 1) + 1, Nt))
    @test_nowarn u[rand(1:Ny), rand(1:((Nz >> 1) + 1)), rand(1:Nt)]
    @test_nowarn u[rand(1:Ny), rand(1:((Nz >> 1) + 1)), rand(1:Nt)] = 2.0
end

@testset "Scalar field broadcasting     " begin
    # construct field
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
    u = RPCFScalarField(g, rand(ComplexF64, Ny, (Nz >> 1) + 1, Nt))

    @test u .+ Complex(1.0, 2.0) == u + RPCFScalarField(g, Complex(1.0, 2.0)*ones(Ny, (Nz >> 1) + 1, Nt))
end

@testset "Scalar field transform        " begin

end

@testset "Scalar field multiplication   " begin

end

@testset "Scalar field vector calculus  " begin

end

@testset "Scalar field inner product    " begin
    # initialise grid variables
    Ny = 64; Nz = 64; Nt = 64
    y = chebpts(Ny)
    Dy = chebdiff(Ny)
    Dy2 = chebddiff(Ny)
    ws = chebws(Ny)
    ω = abs(rand())
    β = abs(rand())

    # initialise grid
    grid = Grid(y, Nz, Nt, Dy, Dy2, ws, ω, typeof(ω)(β))

    # initialise function
    func1(y, z, t) = (1 - y^2)*exp(cos(β*z))*cos(sin(ω*t))
    func2(y, z, t) = cos(π*y)*(1 - y^2)*exp(sin(β*z))*(cos(ω*t)^2)
end
