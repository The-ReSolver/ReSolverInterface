# This file contains tests for a basic concrete implementation of
# AbstractScalarField.

# TODO: add the other required methods
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

@testset "Scalar field construction" begin

end

@testset "Scalar field array interface" begin

end

@testset "Scalar field broadcasting" begin

end

@testset "Scalar field grid methods" begin

end

@testset "Scalar field differentiation" begin

end

@testset "Scalar field inner-product and norm" begin

end
