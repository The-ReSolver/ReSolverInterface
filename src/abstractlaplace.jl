# This file contains the abstract definitions required to define the behaviour
# of a Laplace operator, and how it can be used to solve Poisson problems.

"""
    AbstactLaplace

A Laplace operator over a domain.
"""
abstract type AbstractLaplace end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constructor method
(::Type{<:AbstractLaplace})(::AbstractGrid) = throw(NotImplementedError())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# solution methods
"""
    solve!(ϕ::AbstractArray, Δ::AbstractLaplace, f::AbstractArray) -> AbstractArray

TODO: explanation here!
"""
solve!(phi::AbstractArray, Δ::AbstractLaplace, rhs::AbstractArray) = ldiv!(phi, Δ, rhs)

"""
    solve!()
"""
solve!(::AbstractArray, ::AbstractLaplace, ::AbstractArray, ::NTuple{N, AbstractArray}) where {N} = throw(NotImplementedError())
