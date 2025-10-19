"""
    Statistical properties of the escape rate EWS in a stationary, coarse parameter sweep in a quadratic saddle-node.

Characterise the variance of the escape rate early-warning signal computed from the potential reconstruction by running multiple, stationary ensemble simulations at differe fixed parameter values.
As the bifurcation parameter approaches its critical value the authentic early-warning (the one computed analytically from the ground truth potential U(x)) increases and the reconstructions becomes better which is signalled by an ensemble mean closer to the real value and an ensemble variance drastically decreasing.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the ensemble's sample paths
        printstyled("Computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:length(ensemble.state)
                # Check for tipping
                tipping = find_tipping(ensemble.state[n], check = 0.010, criterion = 0.010, verbose=false)
                tipping.check && (global Nx += 1)

                # Extract the trajectory until the tipping and center it 
                t = ensemble.time
                u = (ensemble.state[n])[1:tipping.index] .- sqrt(-μ)
                tipping.check && display(t[tipping.index])

                # Solve the nonlinear least-squares problem to fit a local cubic potential
                solution = fit_potential(u, n_bins=Nb, noise=σ, optimiser=β, attempts=Na)
                push!(solutions, solution.fit)

                # Perform postprocessing analysis on the solutions
                Vs = analyse(solution.fit)

                #=
                # Plot and export the potential reconstruction 
                include("./scripts/figs.jl")
                plot_solutions(t, u, Vs, n)
                =#
        end

        # Export the csv file
        export_data()
end

# Execute the main
#main()
# Plot the ews
plot_ews()
