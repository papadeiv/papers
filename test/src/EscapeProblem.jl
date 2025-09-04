module EscapeProblem 

# Import packages
using LinearAlgebra, Polynomials, Integrals
using StatsBase, LsqFit

# Import utility functions
include("../utils/reconstruction_quasipotential.jl")
include("../utils/large_deviation_principle.jl")

# Export namespaces
export fit_distribution, get_normalisation_constant, fit_potential, shift_potential
export estimate_escape 

end # module
