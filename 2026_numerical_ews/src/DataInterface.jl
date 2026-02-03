module DataInterface 

# Import packages
using Tables, CSV, MAT, DataFrames
using DocStringExtensions

# Import utility functions
include("../utils/input.jl")
include("../utils/output.jl")

# Export namespaces
export readin, writeout

end # module
