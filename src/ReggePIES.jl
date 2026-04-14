module ReggePIES

using StaticArrays
using Graphs
using LinearAlgebra

include("combinatorics.jl")
include("state.jl")
include("parallel.jl")
include("solver.jl")
include("residuals.jl")
include("evolver.jl")
include("state.jl")
include("localmove.jl")
include("geometry.jl")
include("gauge.jl")
include("diagnostics.jl")
include("initdata.jl")
include("defects.jl")



end


