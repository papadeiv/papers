# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Iterate over the parameters grid 
        for (parameter_index, μ) in enumerate(μ_set)
                # Solve the differential equation 
                solution = evolve(f, η, Λ, [x0, μ], stepsize=dt, steps=Nt)
                t = solution.time
                x = solution.state[1]

                # Plot the timeseries solution
                lines!(ax, t, x, color = parameter_index, colormap = (:phase,1.0), colorrange = (1,length(μ_set)), linewidth = 3.0)
        end

        # Export the figure
        savefig("bifurcation_timeseries.pdf", fig)
end

# Execute the main
main()
