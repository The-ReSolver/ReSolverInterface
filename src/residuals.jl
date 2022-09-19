# This file contains the fallback methods for the calculation of the residual
# quantities of a flow field such that it can be optimised over.

# ! these methods require rework in main file such that the expected functions
# ! ℜ and dℜ can be computed without intermediate steps !

ℜ(::AbstractVector{<:AbstractScalarField}, otherargs...) = throw(NotImplementedError())

dℜ(::AbstractVectir{<:AbstractScalarField}, otherargs...) = throw(NotImplementedError())

function ℜdℜ(U::AbstractVector{<:AbstractScalarField}, otherargs...)
    ℜ = ℜ(U, otherargs...)
    dℜ = dℜ(U, otherargs...)
    return ℜ, dℜ
end
