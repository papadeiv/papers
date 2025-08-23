# Import packages
using LinearAlgebra, Polynomials, Integrals
using StatsBase, LsqFit

# Include functions
include("../src/fit_distribution.jl")
include("../src/fit_potential.jl")
include("../src/estimate_error.jl")
include("../src/approximate_potential.jl")
include("../src/get_normalisation_constant.jl")
include("../src/invert_equilibrium_distribution.jl")
include("../src/get_stationary_points.jl")
include("../src/gaussian.jl")
