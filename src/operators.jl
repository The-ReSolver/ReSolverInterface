# Definitions for the general operators for the Navier-Stokes equations and the
# variational dynamics.

# FIXME: construct concrete vectorfield so user cannot have their own (maybe I should just remove the options?)
# TODO: remove abstract vector field
struct NavierStokesOperator{N, S}
    cache::Vector{VectorField{N, S}}
    Re_recip::Float64

    NavierStokesOperator(::Type{S}, g::AbstractGrid, Re::Float64) where {S<:AbstractScalarField} = new{length(cache[1]), S}([VectorField(S, g) for _ in 1:2], 1/Re)
end

function (f::NavierStokesOperator{N, S})(N_u::VectorField{N, S}, u::VectorField{N, S}) where {N, S<:AbstractScalarField}
    # assign aliases
    u∇u = f.cache[1]
    Δu  = f.cache[2]

    # compute operator
    @. N_u = -convection!(u∇u, u) + f.Re_recip*laplacian!(Δu, u)

    return N_u
end


struct GradientOperator{N, S}
    cache::Vector{VectorField{N, S}}
    Re_recip::Float64

    GradientOperator(::Type{S}, g::AbstractGrid, Re::Float64) where {S<:AbstractScalarField} = new{length(cache[1]), S}([VectorField(S, g) for _ in 1:3], 1/Re)
end

function (f::GradientOperator{N, S})(M_ur::VectorField{N, S}, u::VectorField{N, S}, r::VectorField{N, S}) where {N, S<:AbstractScalarField}
    # assign aliases
    u∇r = f.cache[1]
    ∇ur = f.cache[2]
    Δr  = f.cache[3]

    # compute operator
    @. M_ur = convection!(u∇r) + convection2!(∇ur, u, r) - f.Re_recip*laplacian!(Δr, r)

    return M_ur
end
