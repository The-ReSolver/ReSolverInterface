# This file contains the abstract type and function definitions required for a
# scalar field.

"""
    AbstractScalarField

A scalar field defined over a finite domain.
"""
abstract type AbstractScalarField{N, T<:Number} <: AbstractArray{T, N} end

# ! required !
(::Type{<:AbstractScalarField})(g::AbstractGrid) = throw(NotImplementedError(g))

# ! required !
# (::Type{<:AbstractScalarField})(g::AbstractGrid, func) = throw(NotImplementedError(g, func))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid methods
grid(u::AbstractScalarField) = throw(NotImplementedError(u))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
# ! required !
Base.parent(u::AbstractScalarField) = throw(NotImplementedError(u))

# ! required !
Base.similar(u::AbstractScalarField, ::Type{T}) where {T} = throw(NotImplementedError(u))

# * optional *
Base.size(u::AbstractScalarField) = size(parent(u))

# * optional *
Base.IndexStyle(::Type{<:AbstractScalarField}) = Base.IndexLinear()

# * optional *
Base.copy(u::AbstractScalarField) = (v = similar(u); v .= u; return v)

# * optional *
# TODO: benchmark to see if indexing like this is any faster than just normal indexing
Base.@propagate_inbounds function Base.getindex(u::AbstractScalarField, I...)
    @boundscheck checkbounds(parent(u), I...)
    @inbounds v = parent(u)[I...]
    return v
end

# * optional *
Base.@propagate_inbounds function Base.setindex!(u::AbstractScalarField, v, I...)
    @boundscheck checkbounds(parent(u), I...)
    @inbounds parent(u)[I...] = v
    return v
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pointwise multiplication over the domain
# * optional *
"""
    mult!(
        uv::AbstractScalarField,
        u::AbstractScalarField,
        v::AbstractScalarField
    ) -> AbstractScalarField

Multiply two scalar fields in a pointwise (local) manner and return the result
in `uv`.
"""
mult!(uv::AbstractScalarField, u::AbstractScalarField, v::AbstractScalarField) = throw(NotImplementedError(uv, u, v))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
# ! required !
"""
    LinearAlgebra.dot(
        p::AbstractScalarField,
        q::AbstractScalarField
    ) -> Float64

Compute the inner-product of two scalar fields and return the result.
"""
LinearAlgebra.dot(u::AbstractScalarField, v::AbstractScalarField) = throw(NotImplementedError(u, v))

# *optional *
"""
    LinearAlgebra.norm(p::AbstractScalarField) -> Float64

Compute the norm of a scalar field and return the result.
"""
LinearAlgebra.norm(p::AbstractScalarField) = sqrt(LinearAlgebra.dot(p, p))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derivative methods
# ! required !
"""
    ddt!(
        dudt::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the time derivative of the scalar field u, overwriting dudt with the
result
"""
ddt!(dudt::S, u::S) where {S<:AbstractScalarField} = throw(NotImplementedError(dudt, u))
