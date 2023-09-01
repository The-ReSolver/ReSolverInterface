# This file contains the concrete type definitions to allow the testing of the
# abstract interface defined in the rest of the package.

# grid on a line
struct MyGrid <: AbstractGrid{Float64, 2}
    A::Matrix{Tuple{Float64, Float64}}
end
ReSolverInterface.points(g::MyGrid) = g.A
function ReSolverInterface.points(g::MyGrid, dim::Int)
    pts = zeros(size(points(g), dim))
    for i in axes(points(g), dim)
        pts[i] = points(g)[i][dim]
    end

    return pts
end

struct MyScalarField <: AbstractScalarField{2, Float64}
    grid::MyGrid
    data::Array{Float64, 2}

    MyScalarField(g::MyGrid, ::Type{T}=Float64) where {T} = new(g, rand(T, size(points(g))))
end
Base.parent(u::MyScalarField) = u.data
ReSolverInterface.grid(u::MyScalarField) = u.grid
ReSolverInterface.mult!(uv::MyScalarField, u::MyScalarField, v::MyScalarField) = (@. uv = u*v; return uv)
ReSolverInterface.grad!(∇u::MyScalarField, u::MyScalarField) = (mul!(∇u, DiffMatrix(points(grid(u), 1), 3, 1), u); return ∇u)
ReSolverInterface.laplacian!(Δu::MyScalarField, u::MyScalarField) = (mul!(Δu, DiffMatrix(points(grid(u), 1), 3, 2), u); return Δu)
ReSolverInterface.ddt!(dudt::MyScalarField, u::MyScalarField) = (dudt .= one.(u); return dudt)
ReSolverInterface.dot(u::MyScalarField, v::MyScalarField) = sum(u.*v)
