module PlottingTools

# Import packages
using CairoMakie, Makie.Colors, LaTeXStrings
using DocStringExtensions

# Increase the definition of the figures
CairoMakie.activate!(; px_per_unit = 2)

# Import utility functions
include("../utils/layouts.jl")
include("../utils/palette.jl")

# Export namespaces
export makefig, savefig 

end # module
