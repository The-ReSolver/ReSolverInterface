# Scaling functionality for the inner-product to improve the condition of the
# optimisation. See: Farazmand (2016), An adjoint-based approach for finding
# invariant solutions of Navier–Stokes equations.

abstract type NormScaling end
struct UniformScaling <: NormScaling end
