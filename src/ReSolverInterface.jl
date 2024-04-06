module ReSolverInterface

using LinearAlgebra

export AbstractGrid, points
export AbstractScalarField, grid, mult!, grad!, laplacian!, ddt!, dot, norm
export VectorField, divergence!, cross!

include("notimplementederror.jl")

include("abstractgrid.jl")
include("abstractscalarfield.jl")
include("projectedscalarfield.jl")
include("vectorfield.jl")
include("broadcasting.jl")
include("vectorcalculus.jl")
include("operators.jl")
include("objective.jl")

end
