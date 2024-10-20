# Definition of the cache object that keeps track of all the optimisation
# variables. In addition it will work as a functor to allow it to directly be
# used to compute the residuals.

# TODO: add interface for fourier transforms

struct Objective{G, PERIOD, S, D, P, B, NSE, VDE}
    cache::Vector{VectorField{D, S}}
    projectedCache::Vector{P}
    base::B
    navierStokesOperator::NSE
    gradientOperator::VDE

    function Objective(grid::G, Re::Real, modes, base::B; include_period::Bool=true, navierStokesOperator=NavierStokesOperator(S, grid, Re), gradientOperator=GradientOperator(S, grid, Re)) where {G<:AbstractGrid, B}
        # initialise cache
        cache = [VectorField(grid, numVelComps(G)) for _ in 1:7]
        projectedCache = [ProjectedField(grid, modes) for _ in 1:3]

        new{G, include_period, eltype(cache[1]), numVelComps(G), typeof(projectedCache), B, typeof(navierStokesOperator), typeof(gradientOperator)}(cache, projectedCache, base, navierStokesOperator, gradientOperator)
    end
end

function (f::Objective{G, PERIOD})(a::ProjectedField{G}) where {G, PERIOD}
    # assign aliases
    u    = f.cache[1]
    dudt = f.cache[2]
    N_u  = f.cache[3]
    r    = f.cache[4]
    s    = f.projectedCache[1]

    # expand velocity coefficients into velocity field
    expand!(u, a)

    # include the base profile
    include_base!(u, f.base)

    # compute local residual
    r .= ddt!(dudt, u) .- f.navierStokesOperator(N_u, u)
    expand!(r, project!(s, r))

    return PERIOD ? (norm(r)^2/(2*volume(grid(u))), frequencyGradient(dudt, r)) : norm(r)^2/(2*volume(grid(u)))
end

function (f::Objective{G})(grad::ProjectedField{G}, a::ProjectedField{G}) where {G}
    # assign aliases
    u    = f.cache[1]
    r    = f.cache[4]
    drdt = f.cache[5]
    M_ur = f.cache[6]
    dRdu = f.cache[7]

    # compute residuals
    output = f(a)

    # compute gradient
    dRdu .= .-ddt!(drdt, r) .+ f.gradientOperator(M_ur, u, r)
    project!(grad, dRdu)

    return output
end

frequencyGradient(dudt, r) = (period(grid(r))/(2π))*dot(dudt, r)

function (f::Objective)(R, G, a::ProjectedField)
    if G === nothing
        return f(a)
    else
        return f(G, a)
    end
end

function (f::Objective)(R, G, x::Vector{Float64})
    a = vectorToField!(f.projectedCache[2], x)
    dRda = f.projectedCache[3]
    if G === nothing
        R = f(a)[1]
    else
        R, dRdω = f(dRda, a)
        fieldToVector!(G, dRda, dRdω)
    end
    return R
end
