# Definitions for a concrete implementation of a projected scalar field built
# upon the abstract interface defined for the generic scalar field.

struct ProjectedField{G<:AbstractGrid, N, T, M} <: AbstractArray{T, N}
    modes::M
    grid::G
    field::Array{T, N}

    ProjectedField(modes::M, grid::G, field::Array{T, N}) where {G, M, T, N} = new{G, N, T, M}(modes, grid, field)
end
ProjectedField(grid::G, modes) where {G<:AbstractGrid} = ProjectedField(modes, grid, projectedField(G, modes))
ProjectedField(u::AbstractScalarField, modes) = ProjectedField(grid(u), modes)
ProjectedField(u::VectorField, modes) = ProjectedField(grid(u), modes)

Base.parent(a::ProjectedField) = a.field
Base.size(a::ProjectedField) = size(parent(a))
Base.similar(a::ProjectedField{G}, ::Type{T}) where {G, T} = ProjectedField(a.modes, a.grid, similar(a.field, T))
Base.ndims(::ProjectedField{G, N}) where {G, N} = N
Base.eltype(::ProjectedField{G, N, T}) where {G, N, T} = T
Base.IndexStyle(::Type{<:ProjectedField}) = IndexLinear()
Base.getindex(a::ProjectedField, i::Int) = parent(a)[i]
Base.setindex!(a::ProjectedField, v, i::Int) = (parent(a)[i] = v; return v)

modes(a::ProjectedField) = a.modes
grid(a::ProjectedField) = a.grid

# ! required !
projectedField(::Type{G}, modes) where {G<:AbstractGrid} = throw(NotImplementedError(G, modes))

# ! required !
expand!(u::VectorField{M, S}, a::ProjectedField{N}) where {M, N, S<:AbstractScalarField{N}} = throw(NotImplementedError(u, a))

# ! required !
project!(a::ProjectedField{N}, u::VectorField{M, S}) where {N, M, S<:AbstractScalarField{N}} = throw(NotImplementedError(a, u))
