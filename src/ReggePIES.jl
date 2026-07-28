module ReggePIES

using StaticArrays
using Graphs
using LinearAlgebra

include("combinatorics.jl")
include("state.jl")
include("parallel.jl")
include("localmove.jl")
include("geometry.jl")
include("defects.jl")
include("residuals.jl")
include("gauge.jl")
include("solver.jl")
include("initdata.jl")
include("diagnostics.jl")
include("evolver.jl")

end

