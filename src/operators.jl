# Definitions for the general operators for the Navier-Stokes equations and the
# variational dynamics.

struct NavierStokesOperator{N, S, C, L}
    cache::Vector{VectorField{N, S}}
    Re_recip::Float64
    conv!::C
    lapl!::L

    NavierStokesOperator(::Type{S}, g::AbstractGrid, N::Int, Re::Float64; conv::C=convection!, lapl::L=laplacian!) where {S<:AbstractScalarField, C, L} = new{N, S, C, L}([VectorField(S, g, N) for _ in 1:2], 1/Re, conv, lapl)
end

function (f::NavierStokesOperator{N, S})(N_u::VectorField{N, S}, u::VectorField{N, S}) where {N, S<:AbstractScalarField}
    # assign aliases
    u∇u = f.cache[1]
    Δu  = f.cache[2]

    # compute operator
    @. N_u = -f.conv!(u∇u, u) + f.Re_recip*f.lapl!(Δu, u)

    return N_u
end


struct GradientOperator{N, S, C1, C2, L}
    cache::Vector{VectorField{N, S}}
    Re_recip::Float64
    conv1!::C1
    conv2!::C2
    lapl!::L

    GradientOperator(::Type{S}, g::AbstractGrid, N::Int, Re::Float64; conv1::C1=convection!, conv2::C2=convection2!, lapl::L=laplacian!) where {S<:AbstractScalarField, C1, C2, L} = new{N, S, C1, C2, L}([VectorField(S, g, N) for _ in 1:3], 1/Re, conv1, conv2, lapl)
end

function (f::GradientOperator{N, S})(M_ur::VectorField{N, S}, u::VectorField{N, S}, r::VectorField{N, S}) where {N, S<:AbstractScalarField}
    # assign aliases
    u∇r = f.cache[1]
    ∇ur = f.cache[2]
    Δr  = f.cache[3]

    # compute operator
    @. M_ur = f.conv1!(u∇r) + f.conv2!(∇ur, u, r) - f.Re_recip*f.lapl!(Δr, r)

    return M_ur
end
