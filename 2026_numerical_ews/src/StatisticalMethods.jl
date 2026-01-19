module StatisticalMethods

# Import packages
using LinearAlgebra, StatsBase, LsqFit
using Polynomials, Integrals
using ProgressMeter, DocStringExtensions

# Import utility functions
include("../utils/inference.jl")
include("../utils/preprocess.jl")
include("../utils/approximation.jl")

# Export namespaces
export fit_distribution, normalise, fit_potential
export build_window, detrend, find_tipping

end # module
