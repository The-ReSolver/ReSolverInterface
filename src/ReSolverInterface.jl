module ReSolverInterface

using LinearAlgebra

export AbstractGrid, AbstractScalarField, AbstractLaplace, VectorField

export points
export getβ, getω, getDy, getDy2
export ddt!, ddy!, ddz!, d2dy2!, d2dz2!
export norm
export specialisevectorfieldconstructor

include("notimplementederror.jl")

include("abstractgrid.jl")
include("abstractscalarfield.jl")
include("vectorfield.jl")


end
