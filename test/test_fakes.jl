# This file contains the concrete type definitions to allow the testing of the
# abstract interface defined in the rest of the package.

# grid on a line
struct MyGrid <: AbstractGrid{Float64, 1}
    A::Vector{Float64}
end
ReSolverInterface.points(g::MyGrid) = g.A

struct MyScalarField <: AbstractScalarField{3, Float64}
    grid::MyGrid
    data::Array{Float64, 3}
end



# # ScalarField{S, T}(::Grid) where {S, T} = ScalarField(zeros(S...), T)
# # ScalarField(::Grid, ::Type{T}=Float64) where {T} = ScalarField(zeros(2, 2, 1), T)
# ScalarField(::Grid) = ScalarField{(2, 2, 2), Float64}(zeros(2, 2, 2))
# ScalarField(::Grid, func::Function) = ScalarField{(2, 2, 2), Float64}(func.(ones(2, 2, 2)))
# Base.parent(u::ScalarField) = u.data
# # specialisevectorfieldconstructor(ScalarField)

# println(ScalarField{1, 1} <: AbstractScalarField)
# println(supertype(ScalarField))
# println(subtypes(AbstractScalarField))

# println(typeof(ScalarField(Grid())))
# println(ScalarField(Grid()))

# println(typeof(ScalarField{(2, 2, 2), Float64}(Grid())))
# println(ScalarField{(2, 2, 2), Float64}(Grid()))
