"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
#include("./scripts/proc.jl")
#include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the ensemble's sample paths
        printstyled("Looping over the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:length(ensemble.state)
                # Extract the trajectory until the tipping and center it 
                t = ensemble.time
                u = (ensemble.state[n])[1:(find_tipping(ensemble.state[n], check = 0.010, criterion = 0.010, verbose=false)).index] .- sqrt(-μ)
                display(u)

        end
end

# Execute the main
main()
