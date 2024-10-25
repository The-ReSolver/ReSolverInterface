module ReSolverInterface

using LinearAlgebra

export AbstractGrid, points, volume
export AbstractScalarField, grid, mult!, grad!, laplacian!, convection!, convection2!, ddt!, dot, norm
export NormScaling
export VectorField, divergence!, cross!
export ProjectedField, expand!, project!
export NavierStokesOperator, GradientOperator

include("notimplementederror.jl")

include("abstractgrid.jl")
include("scaling.jl")
include("abstractscalarfield.jl")
include("vectorfield.jl")
include("projectedfield.jl")
include("broadcasting.jl")
include("vectorcalculus.jl")
include("operators.jl")
include("objective.jl")

end
