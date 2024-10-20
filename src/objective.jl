# Definition of the cache object that keeps track of all the optimisation
# variables. In addition it will work as a functor to allow it to directly be
# used to compute the residuals.

# TODO: add frequency gradient
# TODO: add norm scaling
# TODO: add interface for fourier transforms

struct Objective{G, S, D, P, B, NSE, VDE}
    cache::Vector{VectorField{D, S}}
    projectedCache::P
    base::B
    navierStokesOperator::NSE
    gradientOperator::VDE

    function Objective(grid::G, Re::Real, modes, base::B, navierStokesOperator=NavierStokesOperator(S, grid, Re), gradientOperator=GradientOperator(S, grid, Re)) where {G<:AbstractGrid, B}
        # initialise cache
        cache = [VectorField(grid, numVelComps(G)) for _ in 1:7]
        projectedCache = ProjectedField(grid, modes)

        new{G, eltype(cache[1]), numVelComps(G), typeof(projectedCache), B, typeof(navierStokesOperator), typeof(gradientOperator)}(cache, projectedCache, base, navierStokesOperator, gradientOperator)
    end
end

function (f::Objective{G})(a::ProjectedField{G}) where {G}
    # assign aliases
    u    = f.cache[1]
    dudt = f.cache[2]
    N_u  = f.cache[3]
    r    = f.cache[4]
    s    = f.projectedCache

    # expand velocity coefficients into velocity field
    expand!(u, a)

    # include the base profile
    include_base!(u, f.base)

    # compute local residual
    r .= ddt!(dudt, u) .- f.navierStokesOperator(N_u, u)
    expand!(r, project!(s, r))

    return norm(r)^2/(2*volume(grid(u)))
end

function (f::Objective{G})(grad::ProjectedField{G}, a::ProjectedField{G}) where {G}
    # assign aliases
    u    = f.cache[1]
    r    = f.cache[4]
    drdt = f.cache[5]
    M_ur = f.cache[6]
    dRdu = f.cache[7]

    # compute residuals
    R = f(a)

    # compute gradient
    dRdu .= .-ddt!(drdt, r) .+ f.gradientOperator(M_ur, u, r)
    project!(grad, dRdu)

    return R
end

function (f::Objective)(F, G, a::ProjectedField)
    if G === nothing
        return f(a)
    else
        return f(G, a)
    end
end
