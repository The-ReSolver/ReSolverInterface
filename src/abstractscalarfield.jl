# This file contains the abstract type and function definitions required for a
# scalar field.

# * the fast convolution will be achieved as follows:
# *     - concrete implementations of scalar fields will store both the spectral and physical representations of the field
# *     - the transforms between the two will be stored in the object as field (or a field of the grid itself stored as a field)
# *     - this allows the object to be passed and the transform to be performed without initialising any unnecessary arrays

"""
    AbstractScalarField

A scalar field defined over a finite domain.
"""
abstract type AbstractScalarField{N, T<:Number} <: AbstractArray{T, N} end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constructor methods
# ! required !
(::Type{<:AbstractScalarField})(g::AbstractGrid) = throw(NotImplementedError(g))

# * optional *
# TODO: test this, also check no memory assignment is done in the second point
function (u::Type{<:AbstractScalarField})(g::AbstractGrid, f)
    field = u(g)
    field .= f.(points(g))
    return field
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid methods
grid(u::AbstractScalarField) = throw(NotImplementedError(u))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
# ! required !
Base.parent(u::AbstractScalarField) = throw(NotImplementedError(u))

# * optional *
Base.size(u::AbstractScalarField) = size(parent(u))

# * optional *
Base.IndexStyle(::Type{<:AbstractScalarField}) = Base.IndexLinear()

# * optional *
Base.similar(u::AbstractScalarField{S, T}, ::Type{P}=T) where {S, T, P} = typeof(u)(get_grid(u), P)

# * optional *
# TODO: does this work if a different constructor is defined
Base.copy(u::AbstractScalarField) = typeof(u)(copy(parent(u)))

# * optional *
# TODO: benchmark to see if indexing like this is any faster than just normal indexing
Base.@propagate_inbounds function Base.getindex(u::AbstractScalarField, I...)
    @boundscheck checkbounds(parent(u), I...)
    @inbounds v = parent(u)[I...]
    return v
end

# * optional *
Base.Base.@propagate_inbounds function Base.setindex!(u::AbstractScalarField, v, I...)
    @boundscheck checkbounds(parent(u), I...)
    @inbounds parent(u)[I...] = v
    return v
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# broadcasting
# TODO: this block of code needs some serious testing
# * optional *
const AbstractScalarFieldStyle = Broadcast.ArrayStyle{AbstractScalarField}
Base.BroadcastStyle(::Type{<:AbstractScalarField}) = Broadcast.ArrayStyle{<:AbstractScalarField}()
Base.similar(bc::Base.Broadcast.Broadcasted{AbstractScalarFieldStyle}, ::Type{T}) where {T} = similar(find_field(bc), T)

find_field(bc::Base.Broadcast.Broadcasted) = find_field(bc.args)
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args))
find_field(a::AbstractScalarField, rest) = a
find_field(::Any, rest) = find_field(rest)
find_field(x) = x
find_field(::Tuple{}) = nothing


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
# vector calculus methods
# ! required !
"""
    grad!(
        ∇u::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the gradient of the scalar field u, overwriting ∇u with the result.
"""
grad!(∇u::AbstractScalarField{N}, u::AbstractScalarField{N}) where {N} = throw(NotImplementedError(∇u, u))

"""
    laplacian!(
        Δu::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the Laplacian of the scalar field u, overwriting Δu with the result.
"""
laplacian!(Δu::AbstractScalarField{N}, u::AbstractScalarField{N}) where {N} = throw(NotImplementedError(Δu, u))

"""
    ddt!(
        dudt::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the time derivative of the scalar field u, overwriting dudt with the
result
"""
ddt!(dudt::AbstractScalarField{N}, u::AbstractScalarField{N}) where {N} = throw(NotImplementedError(dudt, u))


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
LinearAlgebra.norm(p::AbstractScalarField) = LinearAlgebra.dot(p, p)
