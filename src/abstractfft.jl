# This file contains the definitions to allow a user to easily define an
# implementation of a Fourier transform for their desired field type using
# the FFTW julia interface.

const ESTIMATE = FFTW.ESTIMATE
const EXHAUSTIVE = FFTW.EXHAUSTIVE
const MEASURE = FFTW.MEASURE
const PATIENT = FFTW.PATIENT
const WISDOM_ONLY = FFTW.WISDOM_ONLY
const NO_TIMELIMIT = FFTW.NO_TIMELIMIT

abstract type AbstractFFTPlan end

(::Type{<:AbstractFFTPlan})(::AbstractGrid; flags::UInt32=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) = throw(NotImplementedError())
(f::Type{<:AbstractFFTPlan})(u::AbstractScalarField; flags=EXHAUSTIVE, timelimit=NO_TIMELIMIT) = f(grid(u), flags=flags, timelimit=timelimit)

function (f::AbstractFFTPlan)() end


abstract type AbstractIFFTPlan end

(::Type{<:AbstractIFFTPlan})(::AbstractGrid; flags::UInt32=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) = throw(NotImplementedError())
(f::Type{<:AbstractIFFTPlan})(u::AbstractScalarField; flags=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) = f(grid(u), flags=flags, timelimit=timelimit)

function (f::AbstractIFFTPlan)() end
