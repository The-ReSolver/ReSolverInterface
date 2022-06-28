# This file contains the abstract type and function definitions required for a
# scalar field.

# TODO: add parametric typing stuff to functions!
# TODO: continue filling out the required methods (can I think of any more?)
# TODO: do I make the vector stuff abstract or keep it purely concrete?
# TODO: do I need abtract implementations for the projections and other operators?

"""
    AbstractScalarField

A scalar field defined over a finite domain.
"""
abstract type AbstractScalarField{Ny, Nz, Nt} <: AbstractArray{T, 3} end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# domain/grid methods
"""
    getβ(u::AbstractScalarField) -> Float64

Return the spanwise fundamental wavenumber of the domain.
"""
function getβ(u::AbstractScalarField)
    error("Missing concrete method for getβ")
end

"""
    getDy(u::AbstractScalarField) -> AbstractMatrix

Return the first derivative operator of a given field.
"""
function getDy(u::AbstractScalarField)
    error("Missing concrete method for getDy")
end

"""
    getDy2(u::AbstractScalarField) -> AbstractMatrix

Return the second derivative operator of a given field.
"""
function getDy2(u::AbstractScalarField)
    return getDy(u)*getDy(u)
end

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
function mult!(uv::AbstractScalarField, u::AbstractScalarField, v::AbstractScalarField)
    error("Missing concrete method for mult!")
end


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
function ddt!(dudt::AbstractScalarField, u::AbstractScalarField)
    error("Missing concrete method for ddt!")
end

"""
    ddy!(
        dudy::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the wall-normal derivative of a scalar field `u` and return the result
in `dudy`.
"""
function ddy!(dudy::AbstractScalarField, u::AbstractScalarField)
    error("Missing concrete method for ddy!")
end

"""
    ddz!(
        dudz::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the spanwise derivative of a scalar field `u` and return the result in
`dudx`.
"""
function ddz!(dudz::AbstractScalarField, u::AbstractScalarField)
    error("Missing concrete method for ddz!")
end

"""
    d2dy2!(
        d2udy2::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the second wall-normal derivative of a scalar field `u` and return the
result in `d2udy2`.
"""
function d2dy2!(d2udy2::AbstractScalarField, u::AbstractScalarField)
    ddy!(d2udy2, ddy!(d2udy2, u))
end

"""
    d2dz2!(
        d2udz2::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the second spanwise derivative of a scalar field `u` and return the
result in `d2udz2`.
"""
function d2dz2!(d2udz2::AbstractScalarField, u::AbstractScalarField)
    ddz!(d2udz2, ddz!(d2udz2, u))
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
"""
    LinearAlgebra.dot(
        p::AbstractScalarField,
        q::AbstractScalarField
    ) -> Float64

Compute the inner-product of two scalar fields and return the result.
"""
function LinearAlgebra.dot(u::AbstractScalarField, v::AbstractScalarField)
    error("Missing concrete method for inner-product of scalar fields")
end

"""
    LinearAlgebra.norm(p::AbstractScalarField) -> Float64

Compute the norm of a scalar field and return the result.
"""
LinearAlgebra.norm(p::AbstractScalarField) = LinearAlgebra.dot(p, p)