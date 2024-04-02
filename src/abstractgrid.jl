# This file contains the definitions required to define an abstract grid type.

# TODO: implementation is fucked since grids are represented with tuples of numbers at each point

abstract type AbstractGrid{T<:Real, N} end

"""
    points(g::AbstractGrid) -> AbstractArray{<:Real}

Return an array representing the location of the points making up the grid.
"""
points(g::AbstractGrid) = throw(NotImplementedError(g))

"""
    volume(g::AbstractGrid) -> Float64

Return the volume (or generally metric size) of the domain defined by the grid.
"""
volume(g::AbstractGrid) = throw(NotImplementedError(g))

"""
    ==(g1::AbstractGrid, g2::AbstractGrid) -> Bool

Compare two grids to see if they are the same discrete representation of the
same domain.
"""
Base.:(==)(g1::AbstractGrid, g2::AbstractGrid) = points(g1) == points(g2)
