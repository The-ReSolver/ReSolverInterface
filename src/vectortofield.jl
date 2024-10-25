# Functions to convert a projected field to a vector for optimisation, and
# vice versa.

function vectorToField!(a::ProjectedField, x::Vector{Float64})
    for i in eachindex(a)
        a[i] = x[i] + 1im*x[i + 1]
    end
    setFrequency!(a.grid, x[end])
    return a
end

function fieldToVector!(x::Vector{Float64}, a::ProjectedField)
    for (i, coeff) in enumerate(a)
        x[i]     = real(coeff)
        x[i + 1] = imag(coeff)
    end
    x[end] = frequency(a.grid)
    return x
end
fieldToVector(a::ProjectedField) = fieldToVector!(zeros(2*prod(size(a)) + 1), a)
