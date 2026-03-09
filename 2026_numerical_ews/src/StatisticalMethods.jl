module StatisticalMethods

# Import packages
using LinearAlgebra, StatsBase, LsqFit
using Polynomials, Integrals
using ProgressMeter, DocStringExtensions

# Import Python package for empirical mode decomposition
using PyCall
PyEMD = pyimport("PyEMD")

# Import utility functions
include("../utils/inference.jl")
include("../utils/preprocess.jl")
include("../utils/approximation.jl")

using Tables, CSV, MAT, DataFrames
include("../utils/output.jl")

# Export namespaces
export fit_distribution, normalise, fit_potential
export build_window, detrend, find_tipping

end # module
