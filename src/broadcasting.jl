# This file contains the interface definitions required to get broadcasting
# as desired for concrete scalar and vector fields. Mainly this is to ensure
# that broadcasting on vector fields propagates into the underlying scalar
# field.

Base.BroadcastStyle(u::Type{<:AbstractScalarField}) = Base.Broadcast.ArrayStyle{u}()
Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{S}}, ::Type{T}) where {T, S<:AbstractScalarField} = similar(find_field(bc), T)

Base.BroadcastStyle(q::Type{<:VectorField{N}}) where {N} = Base.Broadcast.ArrayStyle{q}()
Base.similar(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{V}}, ::Type{T}) where {T, N, V<:VectorField{N}} = similar(find_field(bc), T)

Base.BroadcastStyle(::Base.Broadcast.ArrayStyle{V}, ::Base.Broadcast.ArrayStyle{S}) where {V<:VectorField, S<:AbstractScalarField} = Base.Broadcast.ArrayStyle{V}()

find_field(bc::Base.Broadcast.Broadcasted) = find_field(bc.args)
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args))
find_field(a::VectorField, rest) = a
# find_field(a::Base.Broadcast.Extruded{<:VectorField}, rest) = a.x
find_field(a::AbstractScalarField, rest) = a
find_field(::Any, rest) = find_field(rest)
find_field(x) = x
find_field(::Tuple{}) = nothing

function Base.copy(bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{V}}) where {N, V<:VectorField{N}}
    dest = similar(bc, eltype(find_field(bc)[1]))
    for i in 1:N
        copyto!(dest[i], unpack(bc, i))
    end
    return dest
end

function Base.copyto!(dest::VectorField{N}, bc::Base.Broadcast.Broadcasted{Base.Broadcast.ArrayStyle{V}}) where {N, V<:VectorField{N}}
    for i in 1:N
        copyto!(dest[i], unpack(bc, i))
    end
    return dest
end

@inline unpack(bc::Broadcast.Broadcasted, i) = Broadcast.Broadcasted(bc.f, _unpack(i, bc.args))
@inline unpack(x::Any,                    i) = x
@inline unpack(x::VectorField,    i) = x.elements[i]

@inline _unpack(i, args::Tuple) = (unpack(args[1], i), _unpack(i, Base.tail(args))...)
@inline _unpack(i, args::Tuple{Any}) = (unpack(args[1], i),)
@inline _unpack(::Any, args::Tuple{}) = ()
