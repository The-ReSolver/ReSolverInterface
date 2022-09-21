# This file contains the abstract type and function definitions required for a
# scalar field.

# TODO: define matrix interface properly (now needs testing)

# * the fast convolution will be achieved as follows:
# *     - concrete implementations of scalar fields will store both the spectral and physical representations of the field
# *     - the transforms between the two will be stored in the object as field (or a field of the grid itself stored as a field)
# *     - this allows the object to be passed and the transform to be performed without initialising any unnecessary arrays

"""
    AbstractScalarField

A scalar field defined over a finite domain.
"""
abstract type AbstractScalarField{S, T} <: AbstractArray{T, 3} end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constructor methods
# ! required !
(::Type{<:AbstractScalarField})(::AbstractGrid, ::Type{T}) where {T<:Real} = throw(NotImplementedError())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
# ! required !
Base.parent(::AbstractScalarField) = throw(NotImplementedError())

# * optional *
Base.IndexStyle(::Type{<:AbstractScalarField}) = Base.IndexLinear()

# * optional *
Base.similar(u::AbstractScalarField{S, T}, ::Type{P}=T) where {S, T, P} = typeof(U)(get_grid(u), P)

# * optional *
Base.copy(::AbstractScalarField) = (V = similar(U); V .= U; V)

# * optional *
Base.@propagate_inbounds function Base.getindex(u::AbstractArray, I...)
    @boundscheck checkbounds(parent(u), I...)
    @inbounds v = parent(u)[I...]
    return v
end

# * optional *
Base.Base.@propagate_inbounds function Base.setindex!(u::AbstractArray, v, I...)
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
# domain/grid methods
# ! required !
"""
    get_grid(u::AbstractScalarField) -> Type{<:AbstractGrid}

Return the underlying grid on which the given field is defined.
"""
get_grid(u::AbstractScalarField) = throw(NotImplementedError())

# * optional *
"""
    getβ(u::AbstractScalarField) -> Float64

Return the spanwise fundamental wavenumber of a given field.
"""
getβ(u::AbstractScalarField) = getβ(get_grid(u))

# * optional *
"""
    getω(u::AbstractScalarField) -> Float64

Return the fundamental frequency of a given field.
"""
getω(u::AbstractScalarField) = getω(get_grid(u))

# * optional *
"""
    getDy(u::AbstractScalarField) -> AbstractMatrix

Return the first derivative operator of a given field.
"""
getDy(u::AbstractScalarField) = getDy(get_grid(u))

# * optional *
"""
    getDy2(u::AbstractScalarField) -> AbstractMatrix

Return the second derivative operator of a given field.
"""
getDy2(u::AbstractScalarField) = getDy(u)*getDy(u)


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
mult!(::AbstractScalarField, ::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derivative methods
# ! required !
"""
    ddt!(
        dudt::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the time derivative of a scalar field `u` and return the result in
`dudt`.
"""
ddt!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

# ! required !
"""
    ddy!(
        dudy::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the wall-normal derivative of a scalar field `u` and return the result
in `dudy`.
"""
ddy!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

# ! required !
"""
    ddz!(
        dudz::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the spanwise derivative of a scalar field `u` and return the result in
`dudx`.
"""
ddz!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

# * optional *
"""
    d2dy2!(
        d2udy2::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the second wall-normal derivative of a scalar field `u` and return the
result in `d2udy2`.
"""
d2dy2!(d2udy2::AbstractScalarField, u::AbstractScalarField) = ddy!(d2udy2, ddy!(d2udy2, u))

# * optional *
"""
    d2dz2!(
        d2udz2::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the second spanwise derivative of a scalar field `u` and return the
result in `d2udz2`.
"""
d2dz2!(d2udz2::AbstractScalarField, u::AbstractScalarField) = ddz!(d2udz2, ddz!(d2udz2, u))


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
LinearAlgebra.dot(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

# ? leave ?
"""
    LinearAlgebra.norm(p::AbstractScalarField) -> Float64

Compute the norm of a scalar field and return the result.
"""
LinearAlgebra.norm(p::AbstractScalarField) = LinearAlgebra.dot(p, p)
