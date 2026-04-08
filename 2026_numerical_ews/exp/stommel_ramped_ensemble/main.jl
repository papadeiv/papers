"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Loop over the initial parameter values
        for μ0 ∈ M0
                # Compute the initial condition
                equilibria = get_equilibria(f1, f2, μ0)
                z0 = [equilibria.stable[2]; μ0]

                # Solve the ensemble slow-fast SDEs
                ensemble = evolve(f, η, Λ, z0, steps=Nt, stepsize=δt, particles=Ne)

                # Compute the drift of the quasi-steady equilibrium
                qse = [(get_equilibria(f1, f2, μ)).stable[2] for μ in ensemble.parameter]

                # Loop over the ensemble's sample paths
                printstyled("μ=$(μ0): solving the weighted least-squares problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the x- and y-components of the solution, detrend them and compute the observable 
                        x = (ensemble.state[n])[:, 1] .- getindex.(qse, 1)
                        y = (ensemble.state[n])[:, 2] .- getindex.(qse, 2)
                        ψ = x - y 

                        # Check for tipping
                        tipping = find_tipping(ψ, check = 0.010, criterion = 0.050, verbose=false)
                        if tipping.check
                                # Trajectory has tipped => discard it
                        else
                                # Solve the nonlinear least-squares problem
                                solution = fit_potential(ψ, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                                push!(solutions, [ensemble.parameter[end]; solution.fit])

                                # Perform postprocessing analysis on the solutions and plot the results
                                analyse(solution.fit)
                        end
                end

                # Plot the reconstruction of the shifted potential
                include("./scripts/figs.jl")
                plot_solution(ensemble.parameter[end])

                # Export the postprocessed data and vacate the analysis arrays for the next iteration
                export_data()
                empty!(results)
                empty!(solutions)

                # Export the figure
                savefig("stommel/$glb_idx.png", fig1)

                # Update the global index
                println()
                global glb_idx = glb_idx + 1
        end

        # Plot and export the early-warning signal figure
        plot_ews()
end

# Execute the main
main()
