# This file contains the definitions required to perform vector calculus
# operatiosn of the abstract and concrete fields defined elsewhere.

# * optional *
"""
    grad!(
        ∇u::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the gradient of the scalar field u, overwriting ∇u with the result.
"""
grad!(∇u::AbstractVectorField, u::AbstractScalarField) = throw(NotImplementedError(∇u, u))


# ! requried !
divergence!(div_u::AbstractVectorField, u::AbstractVectorField) = throw(NotImplementedError(div_u, u))


# ! required !
"""
    laplacian!(
        Δu::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the Laplacian of the scalar field u, overwriting Δu with the result.
"""
laplacian!(Δu::AbstractScalarField, u::AbstractScalarField) = throw(NotImplementedError(Δu, u))

# ! requried !
laplacian!(Δu::AbstractVectorField, u::AbstractVectorField) = throw(NotImplementedError(Δu, u))


# ! required !
"""
    ddt!(
        dudt::AbstractScalarField,
        u::AbstractScalarField
    ) -> AbstractScalarField

Compute the time derivative of the scalar field u, overwriting dudt with the
result
"""
ddt!(dudt::AbstractScalarField, u::AbstractScalarField) = throw(NotImplementedError(dudt, u))
ddt!(dudt::VectorField, u::VectorField) = ddt!.(dudt, u)
