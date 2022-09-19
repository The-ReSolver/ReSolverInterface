module ReSolverInterface

using LinearAlgebra

export AbstractGrid, AbstractScalarField, AbstractLaplace

export points
export getβ, getω, getDy, getDy2
export ddt!, ddy!, ddz!, d2dy2!, d2dz2!
export norm

include("notimplementederror.jl")
include("abstractgrid.jl")
include("abstractscalarfield.jl")
include("abstractlaplace.jl")
include("residuals.jl")

end
