module SystemAnalysis

# Import packages
using LinearAlgebra, DifferentialEquations
using Roots, ForwardDiff
using DocStringExtensions

# Import utility functions
include("../utils/evolution.jl")
include("../utils/equilibria.jl")

# Export namespaces
export evolve 
export get_equilibria

end # module
