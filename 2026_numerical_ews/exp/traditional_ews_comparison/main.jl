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
        # Build and plot the bifurcation diagram
        equilibria = get_equilibria(f, c0, domain=[-10,10])
        bif_diag = compute_bif_diag(equilibria.stable)
        lines!(ax1, bif_diag[:,1], bif_diag[:,2], color = ifelse.(bif_diag[:,3] .≤ 0, :blue, :red), linewidth = 3.0)

        # Solve the slow-fast SDE
        x0 = [equilibria.stable[1], c0]
        ensemble = evolve(f, η, Λ, x0, endparameter = cf, stepsize = δt, particles = Ne)
        t = ensemble.time
        μ = ensemble.parameter 
        display(length(t))

        # Loop over the ensemble's sample paths
        printstyled("Solving the likelihood estimation problems\n"; bold=true, underline=true, color=:light_magenta)
        @showprogress for (solution_index, x) in enumerate(ensemble.state)
                # Identify the tipping and truncate the solution up to the that point
                tipping_index = find_tipping(x; width = 300, threshold = 2.75, verbose=false) - 100::Integer
                μt = μ[1:tipping_index]
                xt = x[1:tipping_index]

                # Compute the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(xt; alg = "emd", n_modes = 1).residuals

                # Convert the sliding window into an ensemble of suberies
                subseries = preprocess_solution(μt, xt, window_size, verbose = false)

                # Loop over the window strides 
                ews = Matrix{Float64}(undef, length(subseries.trajectories), 5)
                for (n, (parameter, trajectory)) in enumerate(zip(subseries.timesteps, subseries.trajectories))
                        # Extract the parameter value at the end of the window
                        ews[n,1] = parameter[end]

                        # Compute the traditional early-warning signals
                        ews[n,2] = var(trajectory)
                        ews[n,3] = cov(trajectory[1:end-1], trajectory[2:end])/var(trajectory[1:end-1])
                        ews[n,4] = skewness(trajectory)

                        # Compute the escape early-warning signal
                        c = solve_EM_MLE(trajectory)
                        ews[n,5] = compute_ews(c)
                end

                # Plot the results
                lines!(ax1, μ, x, color = (:black, 0.05), linewidth = 1.0)
                #scatter!(ax1, μt[end], x[tipping_index], markersize = 20, color = :yellow, strokewidth = 1.0, strokecolor = :black)
                #lines!(ax1, [μt[end], μt[end]], [minimum(x), maximum(x)], color = (:black,1.5), linestyle = :dash, linewidth = 3.0)
                #lines!(ax2, [μt[end], μt[end]], [minimum(residuals), maximum(residuals)], color = (:black,0.5), linestyle = :dash, linewidth = 3.0)
                lines!(ax2, μt, residuals, color = (:black, 0.05), linewidth = 1.0)
                lines!(ax3L, ews[:,1], ews[:,2], color = (:black, 0.15), linewidth = 1.0)
                lines!(ax4L, ews[:,1], ews[:,3], color = (:black, 0.15), linewidth = 1.0)
                lines!(ax3R, ews[:,1], ews[:,4], color = (:black, 0.15), linewidth = 1.0)
                lines!(ax4R, ews[:,1], ews[:,5], color = (:blue, 0.10), linewidth = 1.0)
        end

        save("../../res/fig/traditional_ews_comparison.pdf", fig)
end

# Execute the main
main()
