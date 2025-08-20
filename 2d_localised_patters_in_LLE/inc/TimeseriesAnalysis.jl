# Import packages
using Statistics, PyCall
using ProgressMeter

# Include functions
include("../src/get_window_parameters.jl")
include("../src/find_tipping.jl")
include("../src/detrend.jl")
