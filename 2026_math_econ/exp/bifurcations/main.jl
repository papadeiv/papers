# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Iterate over the parameters grid 
        for μ in M0
                # Global plotting index
                global plt_idx = 1::Integer
                include("./scripts/figs.jl")

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
                savefig("(τ=$(μ[1]), γ=$(μ[2]), c_tilde=$(μ[3])).png", fig)
        end
end

# Execute the main
main()
