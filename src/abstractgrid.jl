# This file contains the definitions required to define an abstract grid type.

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
    fieldType(::AbstractGrid) -> S<:AbstractScalarField

Return the type of the scalar field associated with a particular grid.
"""
fieldType(::Type{G}) where {G<:AbstractGrid} = throw(NotImplementedError(G))

"""
    numVelComps(::Type{<:AbstractGrid}) -> Int

Return the number of velocity components associated with a complete velocity
field for a given grid.
"""
numVelComps(::Type{G}) where {G<:AbstractGrid} = throw(NotImplementedError(G))

"""
    period(g::AbstractGrid) -> Float64

Return the period of a given temporal grid.
"""
period(g::AbstractGrid) = throw(NotImplementedError(g))

"""
    ==(g1::AbstractGrid, g2::AbstractGrid) -> Bool

Compare two grids to see if they are the same discrete representation of the
same domain.
"""
Base.:(==)(g1::AbstractGrid, g2::AbstractGrid) = points(g1) == points(g2)
