# Definitions for the general operators for the Navier-Stokes equations and the
# variational dynamics.

struct NavierStokesOperator{S, N, V}
    cache::Vector{V}
    Re::Float64

    function NavierStokesOperator(::Type{S}, g::AbstractGrid) where {S<:AbstractScalarField}

    end
end

function (f::NavierStokesOperator{S, N})(N_u::VectorField{N, S}, u::VectorField{N, S}) where {N, S<:AbstracScalarField}
    # assign aliases
    u∇u = f.cache[1]
    Δu  = f.cache[2]

    # compute operator
    @. N_u = -convection!(u∇u, u) + Re*laplacian!(Δu, u)

    return N_u
end


struct GradientOperator{S, N}
    cache
    Re::Float64

    function GradientOperator(::Type{S}, g::AbstractGrid) where {S<:AbstractScalarField}

    end
end

function (f::GradientOperator{S, N})(M_ur::VectorField{N, S}, u::VectorField{N, S}, r::VectorField{N, S}) where {N, S<:AbstractScalarField}

end
