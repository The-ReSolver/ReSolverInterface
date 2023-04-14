module ReSolverInterface

using LinearAlgebra

export AbstractGrid, AbstractScalarField, AbstractLaplace, VectorField

export points
export getβ, getω, getDy, getDy2
export ddt!, ddy!, ddz!, d2dy2!, d2dz2!
export norm
export specialisevectorfieldconstructor

include("notimplementederror.jl")

include("differentialoperators.jl")
include("abstractgrid.jl")
include("abstractscalarfield.jl")
include("vectorfield.jl")

# TODO: label the methods as "required", "optional", or "leave"

# TODO: could add a general approach to differentiation matrices that work with the abstract fields despite their dimensions if I want


end
