using ReSolverInterface
using Test
using LinearAlgebra

using FDGrids

include("test_fakes.jl")

include("test_notimplementederror.jl")
include("test_abstractgrid.jl")
include("test_abstractscalarfield.jl")
include("test_vectorfield.jl")
