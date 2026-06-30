"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")
include("./scripts/figs.jl")
include("./scripts/SN.jl")
include("./scripts/May.jl")
include("./scripts/Stommel.jl")

# Define the main algorithm
function main()
        # Execute the SN script
        SN()

        # Execute the May script
        May()

        # Execute the Stommel script
        Stommel()
end

# Execute the main
main()
