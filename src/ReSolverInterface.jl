module ReSolverInterface

using LinearAlgebra

export getβ, getω, getDy, getDy2
export ddt!, ddy!, ddz!, d2dy2!, d2dz2!
export norm

include("abstractscalarfield.jl")

end
