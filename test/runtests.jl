using ReSolverInterface
using Test
using LinearAlgebra
using FFTW

using FDGrids

include("fake/fake_grid.jl")
include("fake/fake_field.jl")

include("test_notimplementederror.jl")
include("test_abstractgrid.jl")
include("test_abstractscalarfield.jl")
include("test_vectorfield.jl")
