# This file contains the definitions defining the general behaviour of
# differential operators on a discretised grid. A few useful concrete operators
# are provided for some useful base cases.

# The interface defined here is not required to be used by the user since the
# derivartive functions will be dispatched on AbstractArrays, so this code is
# more for convenience of defining new derivatives with as little overhead as
# possible for the user.

abstract type DifferentialOperator{T, N} <: AbstractArray{T, N} end

# TODO: inspect the FDGrids.jl and ChebUtils.jl repos to see what is required by a differentiation matrix
