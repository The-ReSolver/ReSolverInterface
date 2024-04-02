# Definitions for a concrete implementation of a projected scalar field built
# upon the abstract interface defined for the generic scalar field.

# This type technically breaks what is expected from the abstractscalarfield
# interface. This is not a massive problem since its workings should remain
# internal and therefore does not have expose a faulty implementation to the
# user.

struct ProjectedScalarField{S<:AbstractScalarField, N, T, M} <: AbstractScalarField{N, T}
    modes::M
    field::S
end

# ! required !
expand!(u::S, a::ProjectedScalarField{N, S}) where {N, S<:AbstractScalarField{N}} = throw(NotImplementedError(u, a))
project!(a::ProjectedScalarField{N, S}, u::S) where {N, S<:AbstractScalarField{N}} = throw(NotImplementedError(a, u))

Base.parent(a::ProjectedScalarField) = parent(a.field)
Base.similar(a::ProjectedScalarField, ::Type{T}) where {T} = ProjectedScalarField(a.modes, similar(a.field, T))

# ! required !
LinearAlgebra.dot(u::ProjectedScalarField, v::ProjectedScalarField) = throw(NotImplementedError(u, v))
