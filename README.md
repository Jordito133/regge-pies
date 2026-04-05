# regge-pies
A general-purpose simulator using the Parallelizable Implicit Evolution Scheme (PIES) for Regge calculus.

packages:
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

run:
julia --project=. scripts/run.jl
