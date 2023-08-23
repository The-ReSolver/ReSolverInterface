# This file contains the concrete type definitions to allow the testing of the
# abstract interface defined in the rest of the package.

# grid on a line
struct MyGrid <: AbstractGrid{Float64, 2}
    A::Matrix{Float64}
end
ReSolverInterface.points(g::MyGrid) = g.A

struct MyScalarField <: AbstractScalarField{2, Float64}
    grid::MyGrid
    data::Array{Float64, 2}

    MyScalarField(g::MyGrid) = similar(points(g))
end
ReSolverInterface.grid(u::MyScalarField) = u.grid
ReSolverInterface.mult!(uv::MyScalarField, u::MyScalarField, v::MyScalarField) = (@. uv = u*v; return uv)
ReSolverInterface.grad!(∇u::MyScalarField, u::MyScalarField) = (∇u .= DiffMatrix(points(grid(u)), 3, 1)*u; return ∇u)
ReSolverInterface.laplacian!(Δu::MyScalarField, u::MyScalarField) = (Δu .= DiffMatrix(points(grid(u)), 3, 2)*u; return Δu)
ReSolverInterface.ddt!(dudt::MyScalarField, u::MyScalarField) = (dudt .= one.(u); return dudt)
ReSolverInterface.dot(u::MyScalarField, v::MyScalarField) = sum(u.*v)


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
