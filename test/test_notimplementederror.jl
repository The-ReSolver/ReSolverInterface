@testset "NotImplementedError                   " begin
    # unimplemented method based on abstract type
    f1(a::Number) = throw(ReSolverInterface.NotImplementedError(a))

    # test that abstract method returns error
    @test_throws ReSolverInterface.NotImplementedError f1(2)
    @test_throws "f1(Int64) is missing a concrete implementation" f1(2)

    # concrete method implementation
    f1(::Int) = nothing

    # test that the concrete implementation has the expected behaviour
    @test isnothing(f1(2))
end
