module SystemAnalysis

# Import packages
using DifferentialEquations
using Polynomials, Roots, ForwardDiff
using DocStringExtensions

# Import utility functions
include("../utils/evolve_1d_systems.jl")
include("../utils/evolve_2d_systems.jl")
include("../utils/analyse_system.jl")

# Export namespaces
export evolve_1d, evolve_shifted_1d
export get_equilibria, get_stationary_points

end # module
