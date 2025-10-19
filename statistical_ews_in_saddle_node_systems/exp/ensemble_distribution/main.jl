"""
    Statistical properties of the nonlinear potential reconstruction over an ensemble of i.i.d. stationary particles (i.e. for a fixed parameter value).

We characterise the empirical distribution of the solutions of the nonlinear least-square regression on the cubic saddl-enode potential by assembling the histograms of the polynomial coefficients.
Subsequent quantities related to such coefficients by nonlinear transformations (e.g. equilibria (xs,xu) , potential value at equilibria (V(xs),V(xu)), large deviation principle exp(-ΔV) etc...) are also studied by binning an empirical distribution.
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
                # Extract the trajectory until the tipping and center it 
                t = ensemble.time
                u = (ensemble.state[n])[1:(find_tipping(ensemble.state[n], check = 0.010, criterion = 0.010, verbose=false)).index] .- sqrt(-μ)

                # Solve the nonlinear least-squares problem to fit a local cubic potential
                solution = fit_potential(u, n_bins=Nb, noise=σ, optimiser=β, attempts=Na)
                push!(solutions, solution.fit)

                # Perform postprocessing analysis on the solutions
                Vs = analyse(solution.fit)

                # Only enter this condition for the first 100 particles of the ensemble
                if n < 100
                        # Plot the timeseries and the reconstructed potential 
                        plot_solutions(t, u, Vs, n)
                end
        end

        # Export the solutions figure
        savefig("cubic_1.png", fig1)

        # Plot and export the analysis figures
        plot_results(solutions, results)
        savefig("cubic_2.png", fig3)
        savefig("cubic_3.png", fig6)
        savefig("cubic_4.png", fig14)

        # Export the csv file
        export_data()
end

# Execute the main
main()
