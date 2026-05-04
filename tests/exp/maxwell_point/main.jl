"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the ensemble's sample paths
        printstyled("Looping over the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:convert(Int64, Ne)
                # Extract the trajectory
                t = ensemble.time
                u = ensemble.state[n]

                # Build its empirical distribution
                Nb = convert(Int64, ceil(abs(maximum(u)-minimum(u))/(3.49*std(u)*(length(u))^(-1.0/3.0))))
                bins, pdf = fit_distribution(u, n_bins=Nb+1)

                # Solve the non-linear least-squares problem
                initial_guess = rand(4)
                solution = curve_fit(p, bins, pdf, initial_guess).param

                # Plot and export the result
                include("./scripts/figs.jl")
                plot(bins, pdf, solution, n)
        end
end

# Execute the main
main()
