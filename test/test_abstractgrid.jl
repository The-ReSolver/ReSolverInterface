@testset "Grid construction and methods " begin
    # generate random inputs
    Ny = rand(3:20)
    Nz = rand(3:20)
    Nt = rand(3:20)
    randD = rand([-2, -1, 1, 2])
    y = rand(Float64, Ny)
    D1 = rand(Float32, (Ny, Ny))
    D2 = rand(Float16, (Ny + randD, Ny + randD))
    D_sec = rand(Float32, (Ny, Ny))
    w1 = rand(Int128, Ny)
    w2 = rand(Float16, Ny + randD)
    ω = abs(randn())
    β = abs(rand())

    # test point generation
    g1 = RPCFGrid(y, Nz, Nt, D1, D_sec, w1, ω, β)
    pts = points(g1)
    z = range(0, 2π*(1 - 1/Nz), length = Nz)/β # precision differences in operations
    t = range(0, 2π*(1 - 1/Nt), length = Nt)/ω # mean they aren't exactly equal
    @test collect.(pts) ≈ [[y[ny], z[nz], t[nt]] for ny in 1:Ny, nz in 1:Nz, nt in 1:Nt]

    # test comparison
    g2 = RPCFGrid(y, Nz, Nt + 1, D1, D_sec, w1, ω, β)
    g3 = RPCFGrid(rand(Float64, Ny), Nz, Nt, D1, D_sec, w1, ω, β)
    g4 = RPCFGrid(y, Nz, Nt, rand(Float64, (Ny, Ny)), D_sec, w1, ω, β)
    @test g1 != g2
    @test g1 != g3
    @test g1 == g4
end
