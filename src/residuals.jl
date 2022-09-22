# This file contains the fallback methods for the calculation of the residual
# quantities of a flow field such that it can be optimised over.

# ! these methods require rework in main file such that the expected functions
# ! ℜ and dℜ can be computed without intermediate steps !

ℜ(::AbstractVector{<:AbstractScalarField}, otherargs...) = throw(NotImplementedError())

dℜ(::AbstractVector{<:AbstractScalarField}, otherargs...) = throw(NotImplementedError())
