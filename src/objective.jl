# Definition of the cache object that keeps track of all the optimisation
# variables. In addition it will work as a functor to allow it to directly be
# used to compute the residuals.

struct Objective{S, T, N, D, B, M, NSE, VDE}
    grad::ProjectedField{N, T, M}
    cache::Vector{VectorField{D, S}}
    projectedCache::ProjectedField{N, T, M}
    base::B
    free_mean::Bool
    navierStokesOperator::NSE
    gradientOperator::VDE

    function Objective(::Type{S}, grid::AbstractGrid, D::Int, Re::Real, modes::M, base::B, free_mean::Bool, navierStokesOperator=NavierStokesOperator(S, grid, Re), gradientOperator=GradientOperator(S, grid, Re)) where {S<:AbstractScalarField, M, B}
        # initialise residual gradient output
        grad = ProjectedField(grid, modes)

        # initialise cache
        cache = [VectorField(S, grid, D) for _ in 1:6]
        projectedCache = ProjectedField(grid, modes)

        new{typeof(cache[1][1]), eltype(grad), ndims(grad), D, B, M, typeof(navierStokesOperator), typeof(gradientOperator)}(grad, cache, projectedCache, base, free_mean, navierStokesOperator, gradientOperator)
    end
end

function (f::Objective{S})(a::ProjectedField{S}, compute_grad::Bool=true) where {S}
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
    @. r = ddt!(dudt, u) - f.navierStokesOperator(N_u, u)
    expand!(r, project!(s, r))

    # compute gradient
    if compute_grad
        @. dRdu = ddt!(drdt, r) + f.gradientOperator(M_ur, u, r)
        project!(f.grad, dRdu)
    end

    return f.grad, norm(r)^2/volume(grid(u))
end

function (f::Objective{S})(F, G, a::ProjectedField{S}) where {S}
    G === nothing ? F = f(a, false)[2] : (F = f(a, true)[2]; G .= f.grad)
    return F
end
