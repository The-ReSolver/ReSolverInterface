struct MyField <: AbstractScalarField{2, Float64}
    grid::MyGrid
    data::Matrix{ComplexF64}
end
MyField(g::MyGrid; sizes=(2, 2)) = MyField(g, zeros(ComplexF64, sizes...))
ReSolverInterface.grid(u::MyField) = u.grid
Base.parent(u::MyField) = u.data
Base.similar(u::MyField, ::Type{T}=Float64) where {T} = MyField(grid(u), zero(parent(u)))
ReSolverInterface.mult!(uv::MyField, u::MyField, v::MyField) = (uv .= u.*v; return uv)
ReSolverInterface.grad!(∇u::MyField, u::MyField) = (∇u .= 2 .* u; return ∇u)
ReSolverInterface.laplacian!(Δu::MyField, u::MyField) = (Δu .= u.^2; return Δu)
ReSolverInterface.ddt!(dudt::MyField, u::MyField) = (dudt .= sqrt.(u); return dudt)
LinearAlgebra.dot(u::MyField, v::MyField) = dot(parent(u), parent(v))
