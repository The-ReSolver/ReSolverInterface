module ReSolverInterface

using LinearAlgebra
using FFTW

export AbstractGrid, points
export AbstractScalarField, grid, mult!, grad!, laplacian!, ddt!, dot, norm
export VectorField, specialisevectorfieldconstructor, divergence!, cross!

include("notimplementederror.jl")

include("abstractgrid.jl")
include("abstractscalarfield.jl")
include("vectorfield.jl")
include("broadcasting.jl")
include("vectorcalculus.jl")
include("abstractfft.jl")
include("optcache.jl")


end
