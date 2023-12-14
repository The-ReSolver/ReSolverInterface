# This file contains the definition of the cache object that keeps track of all
# the optimisation variables. In addition it will work as a functor to allow
# it to directly be used to compute the residuals.

# To get a general version of this working I have to allow the user to determine 
# what their Navier-Stokes equation looks like. The cache update functions can
# be kept generic as they rely on the vector calculus operations. To help with
# this a set of methods to extract the cache variables can be created.

abstract type AbstractOptCache end

# NOTE: placeholder so I remember what to do here
dudt(cache::AbstractOptCache) = cache.cache[1]

struct OptCache{S<:AbstractScalarField, A<:AbstractArray, B<:AbstractArray} <:AbstractOptCache
    grad::VectorField{3, S}
    cache::Vector{S}
    base::A
    modes::B
    params::Vector{Float64}
    free_mean::Bool

    # NOTE: constructor can provide default size to cache which can be specialised by the user
end

function (f::OptCache)(a)
    # assign aliases
    # expand velocity field
    # apply base flow

    # update the velocity cache
    _update_velocity!(f)

    # compute the un-projected residual
    @. r_unproj = ns(f)

    # project the residual

    # update the residual cache
    _update_residual!(f)

    # compute the un-projected gradient
    @. grad_unproj = grad(f)

    # project the gradient
    # set the mean to zero if needed
end

function _update_velocity!(cache::AbstractOptCache) end

function _update_residual!(cache::AbstractOptCache) end

ns(::OptCache) = throw(NotImplementedError())
grad(::OptCache) = throw(NotImplementedError())
