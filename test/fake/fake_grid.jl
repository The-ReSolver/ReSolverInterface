struct MyGrid <: AbstractGrid{Float64, 1} end
ReSolverInterface.points(::MyGrid) = [(0, 0, 0)]
ReSolverInterface.fieldType(::Type{<:MyGrid}) = MyField
ReSolverInterface.volume(::MyGrid) = 1.0
ReSolverInterface.numVelComps(::Type{MyGrid}) = 3
