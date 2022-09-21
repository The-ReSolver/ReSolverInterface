# This file contains the concrete implementation of the vector fields based
# based on the abstract scalar field defined elsewhere.

"""
    VectorField{Int, <:AbstractScalarField}([elements])

Subtype of vectors with elements that are subtypes of the AbstractScalarField.
"""
struct VectorField{N<:Integer, S<:AbstractScalarField} <: AbstractVector{S}
    elements::Vector{S}

    # constructor using scalar fields as arguments
    function VectorField(elements::Vararg{<:AbstractScalarField, N}) where {N}
        new{N, typeof(elements[1])}(collect(elements))
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
VectorField(fieldtype::Type{<:AbstractScalarField}, grid::AbstractGrid, N::Int=3) = VectorField([fieldtype(grid) for _ in 1:N]...)

# TODO: a corresponding function one needs to be implemented for the abstract scalar field
VectorField(grid::AbstractGrid, funcs::Vararg{<:Function}) = throw(NotImplementedError())

"""
    specialisevectorfieldconstructor(fieldtype::Type{<:AbstractScalarField})

Utility function to provide a specialised version of the vector field
constructor so that the first argument can be ommitted.

NOTE: only one of these new constructors can be defined for each environment
as otherwise it would lead to method ambiguity in the exact constructor to call
for which field type implied.
"""
function specialisevectorfieldconstructor(fieldtype::Type{<:AbstractScalarField})
    @eval begin
        function VectorField(grid::AbstractGrid, N::Int=3)
            $VectorField($fieldtype, grid, N)
        end
    end
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# misc interface methods
Base.parent(q::VectorField) = q.elements
Base.IndexStyle(::Type{<:VectorField}) = Base.IndexLinear()

Base.getindex(q::VectorField, i::Int) = q.elements[i]
Base.setindex!(q::VectorField, i::Int) = (q.elements[i] = v)

Base.size(::VectorField{N}) where {N} = (N,)
Base.length(::VectorField{N}) where {N} = N

# TODO: check if just broadcasting "similar" over the vector field directly is equivalent
Base.similar(q::VectorField) = VectorField(similar.(q.elements)...)
Base.copy(q::VectorField) = copy.(q)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# broadcasting
# * optional *
const VectorFieldStyle = Broadcast.ArrayStyle{VectorField}
Base.BroadcastStyle(::Type{<:VectorField}) = Broadcast.ArrayStyle{VectorField}()

Base.similar(bc::Base.Broadcast.Broadcasted{VectorFieldStyle}, ::Type{T}) where {T} = similar(find_field(bc), T)

find_field(a::VectorField, rest) = a

@inline function Base.copyto!(dest::VectorField{N}, bc::Base.Broadcast.Broadcasted{VectorFieldStyle}) where {N}
    for i in 1:N
        copyto!(dest.elements[i], unpack(bc, i))
    end

    return dest
end

@inline unpack(bc::Base.Broadcast.Broadcasted, i) = Base.Broadcast.Broadcasted(bc.f, _unpack(i, bc.args))
@inline unpack(x::Any, i) = x
@inline unpack(q::VectorField, i) = q.elements[i]

@inline _unpack(i, args::Tuple) = (unpack(args[1], i), _unpack(i, Base.tail(args))...)
@inline _unpack(i, args::Tuple{Any}) = (unpack(args[1], i),)
@inline _unpack(::Any, args::Tuple{}) = ()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# domain/grid methods
get_grid(q::VectorField) = get_grid(q[1])

getβ(q::VectorField) = getβ(get_grid(q))
getω(q::VectorField) = getω(get_grid(q))
getDy(q::VectorField) = getDy(get_grid(q))
getDy2(q::VectorField) = getDy2(get_grid(q))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derivative methods
ddt!(dudt::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = ddt!.(dudt, u)
ddy!(dudy::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = ddy!.(dudy, u)
ddz!(dudz::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = ddz!.(dudz, u)
d2dy2!(d2udy2::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = d2dy2!.(d2udy2, u)
d2dz2!(d2udz2::VectorField{N, S}, u::VectorField{N, S}) where {N, S} = d2dz2!.(d2udz2, u)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# linear algebra methods
LinearAlgebra.dot(q::VectorField{N}, p::VectorField{N}) where {N} = sum(LinearAlgebra.dot(q[i], p[i]) for i in 1:N)

LinearAlgebra.norm(q::VectorField) = sqrt(LinearAlgebra.dot(q, q))
