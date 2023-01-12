using ReSolverInterface
using Test

# The test suite ensures that a concrete implementation actually implements the
# required methods correctly, and the unit tests ensure (given that concrete
# implementation) that the general methods work as expected.

# ! Put basic concrete implementation here and use test suite to ensure it's valid

struct Grid <: AbstractGrid; A::Matrix{Float64}; end
points(g::Grid) = g.A

struct ScalarField{S, T} <: AbstractScalarField{S, T}; data::AbstractArray{T, 3}; end
# ScalarField{S, T}(::Grid) where {S, T} = ScalarField(zeros(S...), T)
# ScalarField(::Grid, ::Type{T}=Float64) where {T} = ScalarField(zeros(2, 2, 1), T)
ScalarField(::Grid) = ScalarField{(2, 2, 2), Float64}(zeros(2, 2, 2))
ScalarField(::Grid, func::Function) = ScalarField{(2, 2, 2), Float64}(func.(ones(2, 2, 2)))
Base.parent(u::ScalarField) = u.data
# specialisevectorfieldconstructor(ScalarField)

# println(ScalarField{1, 1} <: AbstractScalarField)
# println(supertype(ScalarField))
# println(subtypes(AbstractScalarField))

# println(typeof(ScalarField(Grid())))
# println(ScalarField(Grid()))

# println(typeof(ScalarField{(2, 2, 2), Float64}(Grid())))
# println(ScalarField{(2, 2, 2), Float64}(Grid()))


include("test_notimplementederror.jl")
include("test_abstractgrid.jl")
include("test_abstractscalarfield.jl")
include("test_vectorfield.jl")
