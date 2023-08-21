# This file contains the tests available to the user to make sure their
# concrete implementation follows the rules set out by our abstract interface.

function test_gridinterface(g::AbstractGrid{T, N}) where {T, N}
    # get the points making up the grid
    pts = points(g)

    # check the properties of the points
    pts isa AbstractArray{T, N} || error("Grid points are not correct type!")
end
