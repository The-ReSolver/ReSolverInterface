struct MyGrid <: AbstractGrid{Float64, 1} end
ReSolverInterface.points(g::MyGrid) = [(0, 0, 0)]
