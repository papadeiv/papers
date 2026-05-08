"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
#include("./scripts/plot.jl")

# Define the main algorithm
function main()
        fig = Figure(; size = (900,1200))
        ax1 = Axis(fig[1:2,1:2], limits = (c0, cf, nothing, nothing))
        ax2 = Axis(fig[3:4,1:2], limits = (c0, cf, nothing, nothing))
        ax3L = Axis(fig[5,1])
        ax3R = Axis(fig[5,2])
        ax4L = Axis(fig[6,1])
        ax4R = Axis(fig[6,2])

        # Build the bifurcation diagram of the model
        equilibria = get_equilibria(f, c0, domain=[-10,10])
        bif_diag = compute_bif_diag(equilibria.stable)

        # Solve the slow-fast SDE
        x0 = [equilibria.stable[1], c0]
        solution = evolve(f, η, Λ, x0, endparameter = cf, stepsize = δt)
        t = solution.time
        μ = solution.parameter 
        x = solution.state[1] 

        # Identify the tipping and truncate the solution up to the that point
        tipping_index = find_tipping(x; verbose=true)
        μt = μ[1:tipping_index]
        xt = x[1:tipping_index]

        # Compute the quasi-stationary residuals of the truncated timeseries
        residuals = detrend(xt; alg = "emd", n_modes = 1).residuals

        # Convert the sliding window into an ensemble of suberies
        ensemble = preprocess_solution(μt, xt, window_size)

        # Loop over the ensemble's sample paths
        ews = Matrix{Float64}(undef, length(ensemble.trajectories), 5)
        @showprogress for (n, (parameter, trajectory)) in enumerate(zip(ensemble.timesteps, ensemble.trajectories))
                # Extract the parameter value at the end of the window
                ews[n,1] = parameter[end]

                # Compute the early-warning signals
                ews[n,2] = var(trajectory)
                ews[n,3] = cov(trajectory[1:end-1], trajectory[2:end])/var(trajectory[1:end-1])
                ews[n,4] = skewness(trajectory)
                ews[n,5] = inv(ews[n,3])
        end

        # Plot and export the figure
        lines!(ax1, bif_diag[:,1], bif_diag[:,2], color = ifelse.(bif_diag[:,3] .≤ 0, :blue, :red), linewidth = 3.0)
        lines!(ax1, μ, x, color = :black, linewidth = 2.0)
        lines!(ax1, [μt[end], μt[end]], [minimum(x), maximum(x)], color = (:black,0.5), linestyle = :dash, linewidth = 3.0)
        scatter!(ax1, μt[end], x[tipping_index], markersize = 10, color = :yellow, strokewidth = 1.0, strokecolor = :black)
        lines!(ax2, [μt[end], μt[end]], [minimum(residuals), maximum(residuals)], color = (:black,0.5), linestyle = :dash, linewidth = 3.0)
        lines!(ax2, μt, residuals, color = :black, linewidth = 2.0)
        lines!(ax3L, ews[:,1], ews[:,2], color = :black, linewidth = 1.0)
        lines!(ax4L, ews[:,1], ews[:,3], color = :black, linewidth = 1.0)
        lines!(ax3R, ews[:,1], ews[:,4], color = :black, linewidth = 1.0)
        lines!(ax4R, ews[:,1], ews[:,5], color = :black, linewidth = 1.0)
        save("../../res/fig/traditional_ews_comparison.pdf", fig)
end

# Execute the main
main()
