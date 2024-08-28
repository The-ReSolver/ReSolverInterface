# Definitions for a concrete implementation of a projected scalar field built
# upon the abstract interface defined for the generic scalar field.

struct ProjectedField{N, T, M} <: AbstractArray{T, N}
    modes::M
    field::Array{T, N}
end

# ! required !
ProjectedField(grid::AbstractGrid, modes) = throw(NotImplementedError(grid, modes))

Base.parent(a::ProjectedField) = a.field
Base.size(a::ProjectedField) = size(parent(a))
Base.similar(a::ProjectedField, ::Type{T}) where {T} = ProjectedField(a.modes, similar(a.field, T))
Base.ndims(::ProjectedField{N}) where {N} = N
Base.eltype(::ProjectedField{N, T}) where {N, T} = T
Base.IndexStyle(::Type{<:ProjectedField}) = IndexLinear()
Base.getindex(a::ProjectedField, i::Int) = parent(a)[i]
Base.setindex!(a::ProjectedField, v, i::Int) = (parent(a)[i] = v; return v)

modes(a::ProjectedField) = a.modes

# ! required !
expand!(u::VectorField{M, S}, a::ProjectedField{N}) where {M, N, S<:AbstractScalarField{N}} = throw(NotImplementedError(u, a))

# ! required !
project!(a::ProjectedField{N}, u::VectorField{M, S}) where {N, M, S<:AbstractScalarField{N}} = throw(NotImplementedError(a, u))
