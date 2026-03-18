# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Iterate over the initial conditions
        for u0 in U0
                # Build the discrete-time dynamical system
                galor = DeterministicIteratedMap(f, u0, μ)

                # Extract a forward trajectory
                U, t = trajectory(galor, convert(Int64, N))

                # Plot the 'time-paths'
                plot_paths(U)
        end

        # Export the figure
        savefig("switching.pdf", fig)
end

# Execute the main
main()
