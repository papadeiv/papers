"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        u = ensemble.state[1]
        t = ensemble.time
        μ = ensemble.parameter
        display(length(t))

        tipping_index = find_tipping(u; width = 300, threshold = 2.75, verbose=false) - 100::Integer
        μt = μ[1:tipping_index]
        ut = u[1:tipping_index]

        for n in 1:10
                # Generate a new figure
                include("./scripts/figs.jl")

                # Define the partioned set of values of the regularisation coefficients of the LLS problem
                α_set = LinRange(n-1, n, 11)

                # Plot the sample path
                downsample_index = 40::Integer
                lines!(ax1, μ[1:downsample_index:end], u[1:downsample_index:end], color = :black, linewidth = 1.0)

                # Compute and plot the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(ut; alg = "emd", n_modes = 1).residuals
                lines!(ax2, μt[1:2*downsample_index:end], residuals[1:2*downsample_index:end], color = :black, linewidth = 1.0)

                # Place colorbar for the regularisation coefficient
                Colorbar(fig[4,1:2],
                         limits = (α_set[1], α_set[end]),
                         ticks = α_set,
                         colormap = :viridis,
                         vertical = false,
                         label = L"\alpha",
                         width = Relative(1.0),
                         height = 20
                        )

                # Convert the sliding window across the residuals timeseries into an ensemble of suberies
                subseries = preprocess_solution(μt, residuals, window_size, verbose = false)

                # Loop over the values of the regularisation parameter
                for (reg_index, α) in enumerate(α_set)
                        # Loop over the window strides
                        ews = Matrix{Float64}(undef, length(subseries.trajectories), 3)
                        printstyled("α = $α\n\r"; bold=true, underline=true, color=:light_magenta)
                        @showprogress for (n, (parameter, trajectory)) in enumerate(zip(subseries.timesteps, subseries.trajectories))
                                # Compute the variance and escape early-warning signal
                                signal = compute_ews(trajectory, α=α)
                                ews[n,1] = parameter[end]
                                ews[n,2] = signal.variance 
                                ews[n,3] = signal.escape 
                        end

                        # Plot the EWS 
                        lines!(ax3_L, ews[:,1], ews[:,2], color = (:red, 0.1) , linewidth = 1.0)
                        lines!(ax3_R, ews[:,1], ews[:,3], color = reg_index, colormap = (:viridis, 0.5), colorrange = (1,length(α_set)), linewidth = 1.0)
                end

                # Build and plot the bifurcation diagram
                diagram, bifurcations = compute_bif_diag()
                lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :blue, :red), linewidth = 3.0)
                scatter!(ax1, bifurcations[1,1], bifurcations[1,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)

                # Export the figure
                savefig("RLS/[$(α_set[1]),$(α_set[end])].pdf", fig)
        end
end

# Execute the main
main()
