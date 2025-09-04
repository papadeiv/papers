module StatisticalMethods

# Import packages
using ProgressMeter, StatsBase

# Import utility functions
include("../utils/stationary_processes.jl")
include("../utils/transient_processes.jl")

# Export namespaces
export find_tipping, get_window_parameters, detrend 

end # module
