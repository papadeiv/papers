"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Loop over the parameter values
        printstyled("Looping over the parameter values\n"; bold=true, underline=true, color=:light_blue)
        for (index, μ) in enumerate(μ_set)
                # Compute the initial condition 
                equilibria = get_equilibria(f, μ, domain=[-1,1])
                x0 = [equilibria.stable[2], μ]

                # Redefine the noise level for the Maxwell point
                η(x) = index == 2 ? 1.25*σ : σ

                # Solve the ensemble problem 
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt)
                t = ensemble.time
                u = ensemble.state[1]

                # Compute the sample distribution
                bins, pdf = fit_distribution(u, interval = [-1,1], n_bins=51)

                # Plot the results
                domain = LinRange(-1,1,1000)
                lines!((axes[index])[1], domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                lines!((axes[index])[2], u, t, color = (:red,0.5), linewidth = 1.0)
                barplot!((axes[index])[3], bins, pdf, color = (:red,0.5), strokecolor = :black, strokewidth = 2)
        end

        # Export the figure
        savefig("double_well_potential.pdf", fig)
end

# Execute the main
main()
