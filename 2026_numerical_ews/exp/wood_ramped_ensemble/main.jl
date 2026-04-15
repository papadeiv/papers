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
        for h0 ∈ H0
                #=
                println("H = $(h0): x_stb = $((equilibria.stable[1])[1]), x_uns = $((equilibria.unstable[1])[1])")
                Sn_domain = LinRange(-1.0,1.0,1000)
                St_domain = LinRange(-1.0,1.0,1000)
                SN = [f1(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                ST = [f2(Sn,St,h0) for Sn in Sn_domain, St in St_domain]
                fig = Figure()
                ax = Axis(fig[1,1])
                contour!(ax, Sn_domain, St_domain, SN, levels=[0], color=:red, linewidth=3)
                contour!(ax, Sn_domain, St_domain, ST, levels=[0], color=:blue, linewidth=3)
                save("../../res/fig/$glb_idx.png", fig)
                global glb_idx = glb_idx + 1
                =#

                # Compute the initial condition
                equilibria = get_equilibria(f1, f2, h0, guesses=[[Sn0,St0], [-0.5,0.5], [-0.1,0.6]])
                z0 = [equilibria.stable[1]; h0]

                # Solve the ensemble slow-fast SDEs
                ensemble = evolve(f, η, Λ, z0, steps=Nt, stepsize=δt, particles=Ne)

                # Compute the drift of the quasi-steady equilibrium
                qse = [(get_equilibria(f1, f2, h)).stable[1] for h in ensemble.parameter]

                # Loop over the ensemble's sample paths
                printstyled("μ=$(h0): solving the weighted least-squares problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the x-component of the solution (North-Atlantic salinity) 
                        x = (ensemble.state[n])[:, 1] .- getindex.(qse, 1)

                        # Check for tipping
                        tipping = find_tipping(x, check = 0.010, criterion = 0.050, verbose=false)
                        if tipping.check
                                # Trajectory has tipped => discard it
                        else
                                # Solve the nonlinear least-squares problem
                                solution = fit_potential(x, noise=σ, transformation=[0.0,1.0,8.0])
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
                savefig("wood/$glb_idx.png", fig1)

                # Update the global index
                println()
                global glb_idx = glb_idx + 1
        end

        # Plot and export the early-warning signal figure
        plot_ews()
end

# Execute the main
main()
