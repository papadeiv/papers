module StatisticalMethods

# Import packages
using ProgressMeter, DocStringExtensions
using StatsBase

# Import utility functions
include("../utils/stationary_processes.jl")
include("../utils/transient_processes.jl")

# Export namespaces
export get_window_parameters, detrend 
export find_tipping

end # module
