# This file contains the concrete implementation of the vector fields based
# based on the abstract scalar field defined elsewhere.

"""
    VectorField{Int, <:AbstractScalarField}([elements])

Subtype of vectors with elements that are subtypes of the AbstractScalarField.
"""
struct VectorField{N, S} <: AbstractVector{S}
    elements::Vector{S}

    # constructor using scalar fields as arguments
    function VectorField(elements::Vararg{S, N}) where {S<:AbstractScalarField, N}
        new{N, eltype(elements)}(collect(elements))
    end
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# custom interface methods
# ! required !
include_base!(u::VectorField, base) = throw(NotImplementedError(u, base))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# constructor methods
"""
    VectorField(
            fieldtype::Type{<:AbstractScalarField},
            grid::AbstractGrid,
            N::Int=3
    ) -> VectorField{N, fieldtype}

Construct a vector field using a given type of base scalar field type and a
grid on which it is defined.
"""
VectorField(::Type{F}, grid::AbstractGrid, N::Int=3) where {F<:AbstractScalarField} = VectorField([F(grid) for _ in 1:N]...)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid methods
grid(q::VectorField) = grid(q[1])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
Base.parent(q::VectorField) = q.elements
Base.IndexStyle(::Type{<:VectorField}) = Base.IndexLinear()

Base.getindex(q::VectorField, i::Int) = parent(q)[i]
Base.setindex!(q::VectorField, v, i::Int) = (parent(q)[i] = v; return v)

Base.size(::VectorField{N}) where {N} = (N,)
Base.length(::VectorField{N}) where {N} = N

Base.similar(q::VectorField, ::Type{T}=eltype(q[1])) where {T} = VectorField(similar.(parent(q), T)...)
Base.copy(q::VectorField) = VectorField(copy.(parent(q))...)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vector cross-product
function cross!(v_cross_u::VectorField{3}, v::AbstractVector, u::VectorField{3})
    @. v_cross_u[1] = v[2]*u[3] - v[3]*u[2]
    @. v_cross_u[2] = v[3]*u[1] - v[1]*u[3]
    @. v_cross_u[3] = v[1]*u[2] - v[2]*u[1]

    return v_cross_u
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
LinearAlgebra.dot(q::VectorField{N}, p::VectorField{N}) where {N} = sum(LinearAlgebra.dot(q[i], p[i]) for i in 1:N)
LinearAlgebra.norm(q::VectorField) = sqrt(LinearAlgebra.dot(q, q))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derivative methods
# * optional *
function ddt!(dudt::VectorField{N, S}, u::VectorField{N, S}) where {N, S}
    for i in 1:N
        ddt!(dudt[i], u[i])
    end
    return dudt
end
