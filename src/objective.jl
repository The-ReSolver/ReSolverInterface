# Definition of the cache object that keeps track of all the optimisation
# variables. In addition it will work as a functor to allow it to directly be
# used to compute the residuals.

# TODO: add frequency gradient
# TODO: add norm scaling
# TODO: add interface for fourier transforms

struct Objective{G, S, D, P, B, NSE, VDE}
    grad::P
    projectedCache::P
    cache::Vector{VectorField{D, S}}
    base::B
    free_mean::Bool
    navierStokesOperator::NSE
    gradientOperator::VDE

    function Objective(grid::G, Re::Real, modes, base::B, free_mean::Bool, navierStokesOperator=NavierStokesOperator(S, grid, Re), gradientOperator=GradientOperator(S, grid, Re)) where {G<:AbstractGrid, B}
        # initialise residual gradient output
        grad = ProjectedField(grid, modes)

        # initialise cache
        cache = [VectorField(grid, numVelComps(G)) for _ in 1:7]
        projectedCache = ProjectedField(grid, modes)

        new{G, eltype(cache[1]), numVelComps(G), typeof(grad), B, typeof(navierStokesOperator), typeof(gradientOperator)}(grad, projectedCache, cache, base, free_mean, navierStokesOperator, gradientOperator)
    end
end

function (f::Objective{G})(a::ProjectedField{G}, compute_grad::Bool=true) where {G}
    # assign aliases
    u    = f.cache[1]
    dudt = f.cache[2]
    N_u  = f.cache[3]
    r    = f.cache[4]
    drdt = f.cache[5]
    M_ur = f.cache[6]
    dRdu = f.cache[7]
    s    = f.projectedCache

    # expand velocity coefficients into velocity field
    expand!(u, a)

    # include the base profile
    include_base!(u, f.base)

    # compute local residual
    r .= ddt!(dudt, u) .- f.navierStokesOperator(N_u, u)
    expand!(r, project!(s, r))

    # compute gradient
    if compute_grad
        dRdu .= .-ddt!(drdt, r) .+ f.gradientOperator(M_ur, u, r)
        project!(f.grad, dRdu)
    end

    return f.grad, norm(r)^2/(2*volume(grid(u)))
end

function (f::Objective)(F, G, a::ProjectedField)
    G === nothing ? F = f(a, false)[2] : (F = f(a, true)[2]; G .= f.grad)
    return F
end
