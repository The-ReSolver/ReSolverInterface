# This file contains the definition of the cache object that keeps track of all
# the optimisation variables. In addition it will work as a functor to allow
# it to directly be used to compute the residuals.

# To get a general version of this working I have to allow the user to determine 
# what their Navier-Stokes equation looks like. The cache update functions can
# be kept generic as they rely on the vector calculus operations. To help with
# this a set of methods to extract the cache variables can be created.

struct Objective{T, S<:AbstractScalarField, SP<:AbstractScalarField, B<:AbstractArray, M<:AbstractArray}
    grad::VectorField{3, S}
    cache::Vector{S}
    projectedCache::Vector{SP}
    base::B
    modes::M
    params::Vector{T}
    free_mean::Bool

    function Objective(field::S, modes::M, base::B, params::Vector{T}, free_mean::Bool) where {S<:AbstractScalarField, M, B, T}
        # initialise residual gradient output
        grad = similar(field)

        # initialise cache
        cache = [S(g) for _ in 1:1]
        projectedCache = [S(g, modes) for _ in 1:1]

        new{T, S, typeof(projectedCache), B, M}(grad, cache, projectedCache, base, modes, params, free_mean)
    end
end

function (f::Objective{<:Any, <:Any, SP})(a::SP) where {SP}
    
end

function _update_velocity!(cache) end

function _update_residual!(cache) end

ns(::Objective) = throw(NotImplementedError())
grad(::Objective) = throw(NotImplementedError())
