const ESTIMATE = FFTW.ESTIMATE
const EXHAUSTIVE = FFTW.EXHAUSTIVE
const MEASURE = FFTW.MEASURE
const PATIENT = FFTW.PATIENT
const WISDOM_ONLY = FFTW.WISDOM_ONLY
const NO_TIMELIMIT = FFTW.NO_TIMELIMIT

struct FFTPlan!{Nz, Nt, T}
    plan::FFTW.rFFTWPlan{T, FFTW.FORWARD, false, 3, Tuple{Int, Int}}

    function FFTPlan!(A::Array{T, 3};
                        flags::UInt32=EXHAUSTIVE,
                        timelimit::Real=NO_TIMELIMIT) where {T}
        plan = FFTW.plan_rfft(A, (2, 3); flags=flags, timelimit=timelimit)
        new{size(A, 2), size(A, 3), T}(plan)
    end
end

struct IFFTPlan!{Nz, Nt, T}
    plan::FFTW.rFFTWPlan{Complex{T}, FFTW.BACKWARD, false, 3, Tuple{Int, Int}}
    dummy::Array{Complex{T}, 3}

    function IFFTPlan!(Â::Array{Complex{T}}, Nz;
                        flags=UInt32=EXHAUSTIVE,
                        timelimit::Real=NO_TIMELIMIT) where {T}
        plan = FFTW.plan_brfft(Â, Nz, (2, 3); flags=flags, timelimit=timelimit)
        new{Nz, size(Â, 3), T}(plan)
    end
end





struct RPCFGrid{Ny, Nz, Nt, M<:AbstractMatrix, T} <: AbstractGrid{Float64, 3}
    y::Vector{Float64}
    Dy::M
    Dy2::M
    ws::Vector{Float64}
    β::Float64
    ω::Float64
    FFT!::FFTPlan!{Nz, Nt, T}
    IFFT!::IFFTPlan!{Nz, Nt, T}

    function RPCFGrid(y, Nz, Nt, Dy, Dy2, ws, ω, β, ::Type{T}=Float64; flags::UInt32=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) where {T}
        # error checking
        length(y) == length(ws) || throw(ArgumentError("Quadrature weights not a compatible size!"))
        length(y) == size(Dy, 1) == size(Dy, 2)  || throw(ArgumentError("Differentiation matrix not a compatible size!"))
        length(y) == size(Dy2, 1) == size(Dy2, 2)  || throw(ArgumentError("Differentiation matrix not a compatible size!"))

        # generate Fourier transform plans
        spec = Array{Complex{T}, 3}(undef, length(y), (Nz >> 1) + 1, Nt)
        phys = Array{T, 3}(undef, length(y), Nz, Nt)
        FFT! = FFTPlan!(phys, flags=flags, timelimit=timelimit)
        IFFT! = IFFTPlan!(spec, Nz, flags=flags, timelimit=timelimit)

        # partially initialise object
        new{length(y), Nz, Nt, typeof(Dy), T}(y, Float64.(Dy), Float64.(Dy2), Float64.(ws), Float64(β), Float64(ω), FFT!, IFFT!)
    end
end

function ReSolverInterface.points(g::RPCFGrid{Ny, Nz, Nt}) where {Ny, Nz, Nt}
    # construct array to hold tuple points
    pts = Array{NTuple{3, Float64}, 3}(undef, Ny, Nz, Nt)

    # loop over grid assigning tuple values
    for nt in 1:Nt, nz in 1:Nz, ny in 1:Ny
        pts[ny, nz, nt] = (g.y[ny], (nz - 1)*2π/(g.β*Nz), (nt - 1)*2π/(g.ω*Nt))
    end

    return pts
end






struct RPCFScalarField{Ny, Nz, Nt, T, M} <: AbstractScalarField{3, Float64}
    grid::RPCFGrid{Ny, Nz, Nt, M}
    spec::Array{Complex{T}, 3}
    phys::Array{T, 3}

    function RPCFScalarField(g::RPCFGrid{Ny, Nz, Nt, M}, data::Array{Complex{T}, 3}=zeros(ComplexF64, Ny, (Nz >> 1) + 1, Nt)) where {Ny, Nz, Nt, M, T}
        # generate physical field
        phys = Array{T, 3}(undef, Ny, Nz, Nt)

        new{Ny, Nz, Nt, T, M}(g, data, phys)
    end
end
RPCFScalarField(g::RPCFGrid{Ny, Nz, Nt, M}, ::Type{T}=Float64) where {Ny, Nz, Nt, M, T} = RPCFScalarField(g, zeros(Complex{T}, Ny, (Nz >> 1) + 1, Nt))

ReSolverInterface.grid(u::RPCFScalarField) = u.grid
Base.parent(u::RPCFScalarField, type::Symbol=:spec) = getfield(u, type)
Base.similar(u::RPCFScalarField, ::Type{T}=Float64) where {T} = RPCFScalarField(grid(u), T)




# FFTPlan!(g::RPCFGrid, ::Type{T}=Float64; flags::UInt32=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) where {T} = FFTPlan!(parent(RPCFScalarField(g, T), :phys), flags=flags, timelimit=timelimit)
# FFTPlan!(g::RPCFGrid, ::Type{T}=Float64; flags::UInt32=EXHAUSTIVE, timelimit::Real=NO_TIMELIMIT) where {T} = (u = RPCFScalarField(g, T); IFFTPlan!(parent(u), size(parent(u, :phys), 2), flags=flags, timelimit=timelimit))

# function (f::FFTPlan!{Nz, Nt})(u::RPCFScalarField{<:Any, Nz, Nt}) where {Nz, Nt}
#     # perform transform
#     FFTW.unsafe_execute!(f.plan, parent(u, :phys), parent(û))

#     # normalise
#     parent(u, :phys) .*= (1/(Nz*Nt))

#     return u
# end

# function (f::IFFTPlan!)(u::RPCFScalarField, safe::Bool=true)
#     # perform transform
#     _unsafe_ifft!(f, parent(u, :phys), parent(u), f.dummy, Val(safe))

#     return u
# end

# _unsafe_ifft!(plan, u, û, dummy, ::Val{true}) = FFTW.unsafe_execute!(plan.plan, copyto!(dummy, û), u)
# _unsafe_ifft!(plan, u, û, ::Any, ::Val{false}) = FFTW.unsafe_execute!(plan.plan, û, u)

# function (f::Union{FFTPlan!, IFFTPlan!})(u::VectorField{N}) where {N}
#     for i in 1:N
#         f(u[i])
#     end

#     return u
# end