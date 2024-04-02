# This file contains the definition of the cache object that keeps track of all
# the optimisation variables. In addition it will work as a functor to allow
# it to directly be used to compute the residuals.

# To get a general version of this working I have to allow the user to
# determine what their Navier-Stokes equation looks like. The cache update
# functions can be kept generic as they rely on the vector calculus operations.
# To help with this a set of methods to extract the cache variables can be
# created.

struct Objective{S, T, N, B<:AbstractArray, M<:AbstractArray, NSE, VDE}
    grad::ProjectedScalarField{N, S, T, M}
    cache::Vector{VectorField{3, S}}
    projectedCache::ProjectedScalarField{N, S, T, M}
    base::B
    modes::M
    free_mean::Bool
    navierStokesOperator!::NSE
    gradientOperator!::VDE

    # TODO: if nothing is provided for operators then provide default
    function Objective(::Type{S}, grid::AbstractGrid, modes::M, base::B, navierStokesOperator!, gradientOperator!, free_mean::Bool) where {S<:AbstractScalarField, M, B, T}
        # initialise residual gradient output
        grad = S(g)

        # initialise cache
        cache = [VectorField(S, g) for _ in 1:6]
        projectedCache = ProjectedScalarField(modes, S(g))

        params = convert.(eltype(field), params)

        new{S, eltype(grad), ndims(grad), B, M, typeof(navierStokesOperator!), typeof(gradientOperator!)}(grad, cache, projectedCache, base, modes, params, free_mean, navierStokesOperator!, gradientOperator!)
    end
end

function (f::Objective{S})(a::ProjectedScalarField{S}, compute_grad::Bool=true) where {S}
    # assign aliases
    u    = f.cache[1]
    dudt = f.cache[2]
    N_u  = f.cache[3]
    r    = f.cache[4]
    drdt = f.cache[5]
    M_ur = f.cache[6]
    s    = f.projectedCache[1]

    # expand velocity coefficients into velocity field
    expand!(u, a)

    # include the base profile
    include_base!(u, f.base)

    # compute local residual
    @. r = ddt!(dudt, u) - f.navierStokesOperator!(N_u, u)
    expand!(r, project!(s, r))

    # compute gradient
    if compute_grad
        @. dRdu = ddt!(drdt, r) + f.gradientOperator!(M_ur, u, r)
        project!(f.grad, dRdu)
    end

    return f.grad, norm(s)^2/volume(grid(u))
end
function (f::Objective{S})(F, G, a::ProjectedScalarField{S}) where {S}
    G === nothing ? F = f(a, false)[2] : (F = f(a, true)[2]; G .= f.grad)
    return F
end
# TODO: conversion of field to real vector
(f::Objective)(x::Vector) = f(vectorToVelocityCoefficients!(f.proj_cache[2], x), false)[2]
