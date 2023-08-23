# This file contains the definitions required to define an abstract grid type.

abstract type AbstractGrid{T<:Real, N} end

"""
    points(g::AbstractGrid) -> AbstractArray{<:Real, N}

Return an array representing the location of the points making up the grid.
"""
points(::AbstractGrid) = throw(NotImplementedError())

"""
    ==(g1::AbstractGrid, g2::AbstractGrid) -> Bool

Compare two grids to see if they are the same discrete representation of the
same domain.
"""
Base.:(==)(g1::AbstractGrid, g2::AbstractGrid) = points(g1) == points(g2)
