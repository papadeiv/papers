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
                equilibria = get_equilibria(f, μ0, domain=[-10,10])
                x0 = [equilibria.stable[2], μ0]

                # Solve the ensemble slow-fast SDEs and identify the tipping 
                ensemble = evolve(f, η, Λ, x0, steps=Nt, stepsize=δt, particles=Ne)

                # Compute the drift of the quasi-steady equilibrium
                qse = [(get_equilibria(f, μ, domain=[-10,10])).stable[2] for μ in ensemble.parameter]

                # Loop over the ensemble's sample paths
                print_lock = ReentrantLock()
                counter = Base.Threads.Atomic{Int}(0)
                Threads.@threads for n in 1:length(ensemble.state)
                        # Display the counter status
                        c = atomic_add!(counter, 1)
                        lock(print_lock) do
                                printstyled("\r[μ0 = $μ0]: solving the nonlinear least-squares problems using $(Threads.nthreads()) threads: completion $c/$(length(ensemble.state))"; bold=true, underline=true, color=:light_magenta)
                                flush(stdout)
                        end

                        # Check for tipping
                        tipping = find_tipping(ensemble.state[n], check = 0.050, criterion = 0.050, verbose=false)
                        if tipping.check
                                # Trajectory has tipped => discard it
                        else
                                # Extract and detrend the trajectory 
                                residuals = detrend((ensemble.state[n]), qse = qse).residuals

                                # Solve the nonlinear least-squares problem
                                solution = fit_potential(residuals, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                                push!(solutions, [qse[end]; ensemble.parameter[end]; solution.fit])

                                # Perform postprocessing analysis on the solutions and plot the results
                                analyse(solution.fit, ensemble.parameter[end])
                        end
                end

                # Plot the reconstruction of the shifted potential
                include("./scripts/figs.jl")
                plot_solution(ensemble.parameter[end])

                # Export the postprocessed data and vacate the analysis arrays for the next iteration
                export_data()
                empty!(results)
                empty!(solutions)

                # Update the global index
                println()
                global glb_idx = glb_idx + 1
        end
end

# Execute the main
main()
