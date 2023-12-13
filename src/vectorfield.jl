# This file contains the concrete implementation of the vector fields based
# based on the abstract scalar field defined elsewhere.

"""
    VectorField{Int, <:AbstractScalarField}([elements])

Subtype of vectors with elements that are subtypes of the AbstractScalarField.
"""
struct VectorField{N, S<:AbstractScalarField} <: AbstractVector{S}
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

VectorField(::Type{F}, grid::AbstractGrid, funcs::Vararg{Function}) where {F<:AbstractScalarField} = VectorField([F(grid, funcs[i]) for i in 1:length(funcs)]...)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid methods
grid(q::VectorField) = grid(q[1])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
Base.parent(q::VectorField) = q.elements
Base.IndexStyle(::Type{<:VectorField}) = Base.IndexLinear()

Base.getindex(q::VectorField, i::Int) = q.elements[i]
Base.setindex!(q::VectorField, v, i::Int) = (q.elements[i] = v)

Base.size(::VectorField{N}) where {N} = (N,)
Base.length(::VectorField{N}) where {N} = N

# TODO: base similar method should be defined so there is a complete fallback option
Base.similar(q::VectorField) = similar.(q)
Base.copy(q::VectorField) = copy.(q)

Base.zero(q::VectorField) = zero.(q)
Base.vcat(q::VectorField, p::VectorField) = VectorField(q..., p...)

# TODO: test if broadcasting propogates into the underlying scalar fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# broadcasting
# TODO: separate file for field broadcasting???
# TODO: make broadcasting propogate into the spectral fields depending on input
const VectorFieldStyle = Base.Broadcast.ArrayStyle{VectorField}
Base.BroadcastStyle(::Type{<:VectorField}) = Base.Broadcast.ArrayStyle{VectorField}()
Base.similar(bc::Base.Broadcast.Broadcasted{VectorFieldStyle}, ::Type{T}) where {T} = VectorField(similar.([find_field(bc)...])...)

find_field(a::VectorField, ::Any) = a

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derivative methods
divergence!(div_u::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = throw(NotImplementedError(div_u, u))

laplacian!(Δu::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = throw(NotImplementedError(Δu, u))

ddt!(dudt::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = ddt!.(dudt, u)

function cross!(v_cross_u::VectorField{3, S}, v::AbstractVector, u::VectorField{3, S}) where {S}
    @. v_cross_u[1] = v[2]*u[3] - v[3]*u[2]
    @. v_cross_u[2] = v[3]*u[1] - v[1]*u[3]
    @. v_cross_u[3] = v[1]*u[2] - v[2]*u[1]

    return v_cross_u
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
LinearAlgebra.dot(q::VectorField{N}, p::VectorField{N}) where {N} = sum(LinearAlgebra.dot(q[i], p[i]) for i in 1:N)

LinearAlgebra.norm(q::VectorField) = sqrt(LinearAlgebra.dot(q, q))
