# This file contains a utility to create custom error messages when a method
# has not been implemented.

struct NotImplementedError <: Exception
    method::Symbol
    signature

    NotImplementedError(signature...) = new(stacktrace()[2].func, typeof.(signature))
end

function Base.showerror(io::IO, e::NotImplementedError)
    if length(e.signature) == 1
        print(io, e.method, "(", e.signature[1], ") is missing a concrete implementation!")
    else
        print(io, e.method, e.signature, " is missing a concrete implementation!")
    end
end
