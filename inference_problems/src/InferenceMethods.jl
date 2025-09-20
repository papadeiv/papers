module InferenceMethods

# Import packages
using LinearAlgebra, Polynomials, Integrals
using StatsBase, LsqFit
using ProgressMeter, DocStringExtensions

# Import utility functions
include("../utils/potential_reconstruction.jl")
include("../utils/parameters_estimation.jl")

# Export namespaces
export fit_distribution, fit_potential 

end # module
