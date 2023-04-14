# This file contains the definitions required to define an abstract grid type.

# The grid in its most abstract context only needs to represent what the points
# are that discretise the domain. The information required for the differentiation
# and integration of fields on the grid need to be specified by the user in
# the concrete implementation.

# Or maybe it's best to go even more abstract. Instead of the grid just containing
# the differential operators on a cartesian grid, the differential operators stored can
# be the required vector calculus operations (grad, div, curl, etc.). The exception to
# this will be the time derivative since it is always Fourier transformed for the
# the methods used (assumed stationary flow).

# This would mean a user has to decide how each of the vector calculus operations
# is defined in terms of the individual derivatives. It would also mean the Poisson
# solver could be made concrete since all it requires would be the laplacian from
# the grid, the inhomogeneous data, and boundary data.

abstract type AbstractGrid end

"""
    points(g::AbstractGrid) -> NTuple{3, Vector{Float64}}

Return a tuple of vectors representing the location of the points making up the
grid in each direction.
"""
points(::AbstractGrid) = throw(NotImplementedError())

"""
    ==(g1::AbstractGrid, g2::AbstractGrid) -> Bool

Compare two grids to see if they are the same discrete representation of the
same domain.
"""
Base.:(==)(g1::AbstractGrid, g2::AbstractGrid) = (points(g1) == points(g2))

"""
    Dt(g::AbstractGrid) -> DifferentialOperator

Return the time derivative operator defined on the provided grid.
"""
Dt(::AbstractGrid) = throw(NotImplementedError())

"""
    grad(g::AbstractGrid) -> DifferentialOperator

Return the gradient operator defined on the provided grid.
"""
grad(::AbstractGrid) = throw(NotImplementedError())

"""
    div(g::AbstractGrid) -> DifferentialOperator

Return the divergence operator defined on the provided grid.
"""
div(::AbstractGrid) = throw(NotImplementedError())

"""
    curl(g::AbstractGrid) -> DifferentialOperator

Return the curl operator defined on the provided grid.
"""
curl(::AbstractGrid) = throw(NotImplementedError())

"""
    laplace(g::AbstractGrid) -> DifferentialOperator

Returnt the Laplace operator defined on the provided grid.
"""
laplace(g::AbstractGrid) = dot(div(g), grad(g))
