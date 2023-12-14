# This file contains the concrete implementation of the vector fields based
# based on the abstract scalar field defined elsewhere.

abstract type AbstractVectorField{N, S} <: AbstractVector{S} end

"""
    VectorField{Int, <:AbstractScalarField}([elements])

Subtype of vectors with elements that are subtypes of the AbstractScalarField.
"""
struct VectorField{N, S<:AbstractScalarField} <: AbstractVectorField{N, S}
    elements::Vector{S}

    # constructor using scalar fields as arguments
    function VectorField(elements::Vararg{S, N}) where {S<:AbstractScalarField, N}
        new{N, eltype(elements)}(collect(elements))
    end
end


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

# TODO: see if there is a way to get this working without ambiguities
# VectorField(::Type{F}, grid::AbstractGrid, funcs::Vararg{Function}) where {F<:AbstractScalarField} = VectorField([F(grid, funcs[i]) for i in 1:length(funcs)]...)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid methods
grid(q::AbstractVectorField) = throw(NotImplementedError(q))
grid(q::VectorField) = grid(q[1])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
Base.parent(q::AbstractVectorField) = throw(NotImplementedError(q))
Base.parent(q::VectorField) = q.elements
Base.IndexStyle(::Type{<:AbstractVectorField}) = Base.IndexLinear()

Base.getindex(q::AbstractVectorField, i::Int) = parent(q)[i]
Base.setindex!(q::AbstractVectorField, v, i::Int) = (parent(q)[i] = v; return v)

Base.size(::AbstractVectorField{N}) where {N} = (N,)
Base.length(::AbstractVectorField{N}) where {N} = N

Base.similar(q::AbstractVectorField) = similar.(q)
Base.copy(q::AbstractVectorField) = copy.(q)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# broadcasting
Base.BroadcastStyle(q::Type{<:AbstractVectorField}) = Base.Broadcast.ArrayStyle{q}()
Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{V}}, ::Type{T}) where {T, V<:AbstractVectorField} = VectorField(similar.(find_field(bc).elements)...)

find_field(bc::Base.Broadcast.Broadcasted) = find_field(bc.args)
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args))
find_field(a::AbstractVectorField, ::Any) = a
find_field(a::AbstractScalarField, ::Any) = a
find_field(::Any, rest) = find_field(rest)
find_field(x) = x
find_field(::Tuple{}) = nothing

# @inline function Base.copyto!(dest::VectorField{N}, bc::Base.Broadcast.Broadcasted{VectorFieldStyle}) where {N}
#     for i in 1:N
#         copyto!(dest.elements[i], unpack(bc, i))
#     end

#     return dest
# end

# @inline unpack(bc::Base.Broadcast.Broadcasted, i) = Base.Broadcast.Broadcasted(bc.f, _unpack(i, bc.args))
# @inline unpack(x::Any, i) = x
# @inline unpack(q::VectorField, i) = q.elements[i]

# @inline _unpack(i, args::Tuple) = (unpack(args[1], i), _unpack(i, Base.tail(args))...)
# @inline _unpack(i, args::Tuple{Any}) = (unpack(args[1], i),)
# @inline _unpack(::Any, args::Tuple{}) = ()


function cross!(v_cross_u::AbstractVectorField{3}, v::AbstractVector, u::AbstractVectorField{3})
    @. v_cross_u[1] = v[2]*u[3] - v[3]*u[2]
    @. v_cross_u[2] = v[3]*u[1] - v[1]*u[3]
    @. v_cross_u[3] = v[1]*u[2] - v[2]*u[1]

    return v_cross_u
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
LinearAlgebra.dot(q::AbstractVectorField{N}, p::AbstractVectorField{N}) where {N} = sum(LinearAlgebra.dot(q[i], p[i]) for i in 1:N)

LinearAlgebra.norm(q::AbstractVectorField) = sqrt(LinearAlgebra.dot(q, q))
