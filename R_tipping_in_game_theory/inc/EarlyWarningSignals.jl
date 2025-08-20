# Import packages
using Statistics, StatsBase 

# Define the module
module EarlyWarningSignals

# Include functions
include("./TimeseriesAnalysis.jl")
include("../src/variance.jl")
include("../src/skew.jl")

end
