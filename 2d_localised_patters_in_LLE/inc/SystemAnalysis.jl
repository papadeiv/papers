# Import packages
using DifferentialEquations 
using Roots, ForwardDiff 
using ProgressMeter

# Include functions
include("../src/evolve_forward_1d.jl")
include("../src/evolve_forward_2d.jl")
include("../src/evolve_fixed_1d.jl")
include("../src/evolve_fixed_2d.jl")
include("../src/evolve_ramped_1d.jl")
include("../src/evolve_ramped_2d.jl")
include("../src/evolve_shifted_1d.jl")
include("../src/evolve_ensemble.jl")
include("../src/evolve_ensemble_fixed.jl")
include("../src/get_equilibria.jl")
