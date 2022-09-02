# This file contains the definitions required to define an abstract grid type.

abstract type AbstractGrid end

"""
    points(g::AbstractGrid) -> NTuple{3, Vector{Float64}}

Return a tuple of vectors representing the location of the points making up the
grid in each direction.
"""
points(::AbstractGrid) = throw(NotImplementedError())

"""
    size(g::AbstractGrid) -> NTuple{3, Int}

Return the size of the grid as if it were an array.
"""
Base.size(g::AbstractGrid) = length.(points(g))

"""
    ==(g1::AbstractGrid, g2::AbstractGrid) -> Bool

Compare two grids to see if they are the same discrete representation of the
same domain.
"""
Base.:(==)(g1::AbstractGrid, g2::AbstractGrid) = (points(g1) == points(g2))

"""
    getβ(u::AbstractGrid) -> Float64

Return the spanwise fundamental wavenumber of the domain.
"""
getβ(::AbstractGrid) = throw(NotImplementedError())

"""
    getω(u::AbstractGrid) -> Float64

Return the fundamental frequency of the domain.
"""
getω(::AbstractGrid) = throw(NotImplementedError())

"""
    getDy(u::AbstractGrid) -> AbstractMatrix

Return the first derivative operator of the domain.
"""
getDy(::AbstractGrid) = throw(NotImplementedError())

"""
    getDy2(u::AbstractScalarField) -> AbstractMatrix

Return the second derivative operator of the domain.
"""
getDy2(::AbstractGrid) = throw(NotImplementedError())
