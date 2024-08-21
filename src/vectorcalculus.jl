# This file contains the definitions required to perform vector calculus
# operatiosn of the abstract and concrete fields defined elsewhere. These
# operations are optional, as it is possible to pass your own methods to the
# main operators for more flexible operations. If however, these methods are
# directly extended then no extra input is required.

# * optional *
"""
    grad!(
        ∇u::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the gradient of the scalar field u, overwriting ∇u with the result.
"""
grad!(∇u::V, u::S) where {N, S<:AbstractScalarField, V<:VectorField{N, S}} = throw(NotImplementedError(∇u, u))


# * optional *
"""
    divergence!(
        div_u::AbstractScalarField,
        u::VectorField
    ) -> AbstractScalarField

Compute the divergence of a vector field, overwriting the output into div_u.
"""
divergence!(div_u::S, u::V) where {N, S<:AbstractScalarField, V<:VectorField{N, S}} = throw(NotImplementedError(div_u, u))


# * optional *
"""
    laplacian!(
        Δu::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the Laplacian of the scalar field u, overwriting Δu with the result.
"""
laplacian!(Δu::S, u::S) where {S<:AbstractScalarField} = throw(NotImplementedError(Δu, u))

# * optional *
"""
    laplacian!(
        Δu::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the Laplacian of the vector field u, overwriting Δu with the result.

WARNING: The default implementation assumes cartesian coordinates
"""
function laplacian!(Δu::VectorField{N, S}, u::VectorField{N, S}) where {N, S}
    for i in 1:N
        laplacian!(Δu[i], u[i])
    end
    return Δu
end

# * optional *
"""
    convection!(
        u∇v::VectorField,
        u::VectorField,
        v::VectorField
    ) -> VectorField

Compute the nonlinear convection of a vector field, overwriting the output into
u∇v.
"""
convection!(u∇v::V, u::V, v::V) where {V<:VectorField} = throw(NotImplementedError(u∇v, u, v))
convection!(u∇u::V, u::V) where {V<:VectorField} = convection!(u∇u, u, u)

convection2!(∇uv::V, u::V, v::V) where {V<:VectorField} = throw(NotImplementedError(∇uv, u, v))
