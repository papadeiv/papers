module PlottingTools

# Import packages
using CairoMakie, Makie.Colors
using LaTeXStrings

# Increase the definition of the figures
CairoMakie.activate!(; px_per_unit = 2)

# Import utility functions
include("../utils/figure_layouts.jl")

# Export namespaces
export makefig, savefig 

end # module
