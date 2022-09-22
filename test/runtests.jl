using ReSolverInterface
using Test

# using InteractiveUtils

struct Grid <: AbstractGrid end
struct ScalarField{S, T} <: AbstractScalarField{S, T}; data::AbstractArray{T, 3}; end
# ScalarField{S, T}(::Grid) where {S, T} = ScalarField(zeros(S...), T)
# ScalarField(::Grid, ::Type{T}=Float64) where {T} = ScalarField(zeros(2, 2, 1), T)
ScalarField(::Grid) = ScalarField{(2, 2, 2), Float64}(zeros(2, 2, 2))
ScalarField(::Grid, func::Function) = ScalarField{(2, 2, 2), Float64}(func.(ones(2, 2, 2)))
Base.parent(u::ScalarField) = u.data

# println(typeof(ScalarField(Grid())))
# println(ScalarField(Grid()))

# specialisevectorfieldconstructor(ScalarField)
a = VectorField(ScalarField(Grid()), ScalarField(Grid()))

# println(ScalarField{1, 1} <: AbstractScalarField)
# println(supertype(ScalarField))
# println(subtypes(AbstractScalarField))

# println(typeof(ScalarField{(2, 2, 2), Float64}(Grid())))
# println(ScalarField{(2, 2, 2), Float64}(Grid()))

b = VectorField(ScalarField, Grid(), 3)
c = VectorField(ScalarField, Grid(), x->x^2, x->cos(x))
specialisevectorfieldconstructor(ScalarField)
d = VectorField(Grid(), 4)
e = VectorField(Grid(), x->x^2, x->cos(x))

println(typeof(a))
println(typeof(b))
println(typeof(c))
println(typeof(d))
println(typeof(e))

# TODO: set out the required test sets for the abstract scalar field

@testset "NotImplementedError" begin

end
