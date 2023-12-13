# This file contains the abstract type and function definitions required for a
# scalar field.

# TODO: make this code more generic by removing type information from the higher-order methods

# * the fast convolution will be achieved as follows:
# *     - concrete implementations of scalar fields will store both the spectral and physical representations of the field
# *     - the transforms between the two will be stored in the object as field (or a field of the grid itself stored as a field)
# *     - this allows the object to be passed and the transform to be performed without initialising any unnecessary arrays

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
# broadcasting
# * optional *
const AbstractScalarFieldStyle = Base.Broadcast.ArrayStyle{AbstractScalarField}
Base.BroadcastStyle(u::Type{<:AbstractScalarField}) = Base.Broadcast.ArrayStyle{u}()
Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{S}}, ::Type{T}) where {T<:Real, S<:AbstractScalarField} = similar(find_field(bc), T)
Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{S}}, ::Type{<:Complex{T}}) where {T, S<:AbstractScalarField} = similar(find_field(bc), T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pointwise multiplication over the domain
# ! required !
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
