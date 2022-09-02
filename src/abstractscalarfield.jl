# This file contains the abstract type and function definitions required for a
# scalar field.

# TODO: do I make the vector stuff abstract or keep it purely concrete?
# TODO: do I need abtract implementations for the projections and other operators?
# NOTE: the fast convolution will be achieved as follows:
#           - concrete implementations of scalar fields will store both the spectral and physical representations of the field
#           - the transforms between the two will be stored in the object as field (or a field of the grid itself stored as a field)
#           - this allows the object to be passed and the transform to be performed without initialising any unnecessary arrays 

"""
    AbstractScalarField

A scalar field defined over a finite domain.
"""
abstract type AbstractScalarField{T} <: AbstractArray{T, 3} end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constructor methods
(::Type{<:AbstractScalarField})(::AbstractGrid) = throw(NotImplementedError())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# domain/grid methods
"""
    getβ(u::AbstractScalarField) -> Float64

Return the spanwise fundamental wavenumber of a given field.
"""
getβ(u::AbstractScalarField) = getβ(get_grid(u))

"""
    getω(u::AbstractScalarField) -> Float64

Return the fundamental frequency of a given field.
"""
getω(u::AbstractScalarField) = getω(get_grid(u))

"""
    getDy(u::AbstractScalarField) -> AbstractMatrix

Return the first derivative operator of a given field.
"""
getDy(u::AbstractScalarField) = getDy(get_grid(u))

"""
    getDy2(u::AbstractScalarField) -> AbstractMatrix

Return the second derivative operator of a given field.
"""
getDy2(u::AbstractScalarField) = getDy(u)*getDy(u)

"""
    get_grid(u::AbstractScalarField) -> Type{<:AbstractGrid}

Return the underlying grid on which the given field is defined.
"""
get_grid(u::AbstractScalarField) = u.grid


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pointwise multiplication over the domain
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
"""
    ddt!(
        dudt::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the time derivative of a scalar field `u` and return the result in
`dudt`.
"""
ddt!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

"""
    ddy!(
        dudy::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the wall-normal derivative of a scalar field `u` and return the result
in `dudy`.
"""
ddy!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

"""
    ddz!(
        dudz::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the spanwise derivative of a scalar field `u` and return the result in
`dudx`.
"""
ddz!(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

"""
    d2dy2!(
        d2udy2::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the second wall-normal derivative of a scalar field `u` and return the
result in `d2udy2`.
"""
d2dy2!(d2udy2::AbstractScalarField, u::AbstractScalarField) = ddy!(d2udy2, ddy!(d2udy2, u))

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
"""
    LinearAlgebra.dot(
        p::AbstractScalarField,
        q::AbstractScalarField
    ) -> Float64

Compute the inner-product of two scalar fields and return the result.
"""
LinearAlgebra.dot(::AbstractScalarField, ::AbstractScalarField) = throw(NotImplementedError())

"""
    LinearAlgebra.norm(p::AbstractScalarField) -> Float64

Compute the norm of a scalar field and return the result.
"""
LinearAlgebra.norm(p::AbstractScalarField) = LinearAlgebra.dot(p, p)