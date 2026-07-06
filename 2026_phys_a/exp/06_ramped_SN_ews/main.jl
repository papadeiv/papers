"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        t = ensemble.time
        μ = ensemble.parameter

        # Loop over the ensemble's sample paths
        μ_min, μ_max = -1.0, 0.0
        for (solution_index, solution) in enumerate(ensemble.state)
                # Extract the length of the sample path and truncate the parameter solution to match it
                N_end = length(solution)
                μ_end = μ[1:N_end]

                # Identify the tipping and truncate the solution up to the that point
                tipping_index = find_tipping(solution; width = 300, threshold = 2.75, verbose=false) - 100::Integer
                μt = μ[1:tipping_index]
                xt = solution[1:tipping_index]

                # Compute and plot the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(xt; alg = "emd", n_modes = 1).residuals

                # Convert the sliding window across the residuals timeseries into an ensemble of suberies
                subseries = preprocess_solution(μt, residuals, window_size, verbose = false)

                # Update the maximum value of the parameter in the sliding window
                if (subseries.timesteps[1])[end] > μ_min
                        μ_min = (subseries.timesteps[1])[end]
                end
                # Update the minimum value of the parameter in the sliding window
                if (subseries.timesteps[end])[end] < μ_max
                        μ_max = (subseries.timesteps[end])[end]
                end

                # Plot the sample path and quasi-stationary residuals
                if solution_index ≤ 20
                        downsample_index = 50::Integer
                        lines!(ax1, μ_end[1:downsample_index:end], solution[1:downsample_index:end], color = (:black,0.2), linewidth = 1.0)
                        lines!(ax2, μt[1:2*downsample_index:end], residuals[1:2*downsample_index:end], color = (:black, 0.2), linewidth = 1.0)
                end

                # Loop over the window strides 
                ews = Matrix{Float64}(undef, length(subseries.trajectories), 5)
                printstyled("Progressing trajectory $solution_index of $(convert(Integer, Ne))\r"; bold=true, underline=true, color=:light_magenta)
                for (n, (parameter, trajectory)) in enumerate(zip(subseries.timesteps, subseries.trajectories))
                        # Solve the LLS problem
                        θ = solve_lls(trajectory)

                        # Compute the variance EWS and update the storage matrix
                        variance = var(trajectory)
                        ews[n,1] = parameter[end]
                        ews[n,2] = variance 
                        ews[n,3:5] = θ 
                end

                # Export the EWS timeseries
                writeout(ews, "ews/$solution_index.csv")
        end

        # Compute the number of steps in the filtered subseries
        ews = readin("ews/1.csv")
        mask = (ews[:, 1] .>= μ_min) .& (ews[:, 1] .<= μ_max)
        filtered_ews = ews[mask, :]
        Nt_filtered = size(filtered_ews, 1)
        μ_filtered = filtered_ews[:,1]

        # Loop over the trajectories
        ensemble_θ1 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_θ2 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_θ3 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_variance = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        printstyled("\nComputing ensemble means\n"; bold=true, underline=true, color=:light_magenta)
        @showprogress for n in 1:convert(Integer, Ne)
                # Import the matrix from file
                ews = readin("ews/$n.csv")

                # Filter out the subseries bounded between μ_min and μ_max
                mask = (ews[:, 1] .>= μ_min) .& (ews[:, 1] .<= μ_max)
                filtered_ews = ews[mask, :]

                # Loop over the rows of the filtered matrix
                for (row_index, row) in enumerate(eachrow(filtered_ews))
                        # Update the storage structures
                        ensemble_variance[n,row_index] = row[2]
                        ensemble_θ1[n,row_index] = row[3]
                        ensemble_θ2[n,row_index] = row[4]
                        ensemble_θ3[n,row_index] = row[5]
                end
        end

        # Compute the ensemble means
        ensemble_mean_θ = hcat(vec(mean(ensemble_θ1, dims=1)), vec(mean(ensemble_θ2, dims=1)), vec(mean(ensemble_θ3, dims=1)))
        ensemble_mean_variance = vec(mean(ensemble_variance, dims=1))

        # Compute and plot the ensemble averaged EWS
        ensemble_mean_ews = [compute_ews(θ_mean) for θ_mean in eachrow(ensemble_mean_θ)]
        lines!(ax3_L, μ_filtered, ensemble_mean_variance, color = (:red,0.75), linewidth = 3.0)
        lines!(ax3_R, μ_filtered, ensemble_mean_ews, color = (:red,0.75), linewidth = 3.0)

        # Plot the tipping point and format the axes
        lines!(ax1, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        lines!(ax2, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        lines!(ax3_L, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        lines!(ax3_R, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        ax3_L.limits = (μ_min, 0, 0, 0.02)
        ax3_L.xticks = [μ_min, 0]
        ax3_L.yticks = [0, 0.02]
        ax3_R.limits = (μ_min, 0, 0.25, 1)
        ax3_R.xticks = [μ_min, 0]
        ax3_R.yticks = [0.25,1]

        # Build and plot the bifurcation diagram
        diagram, bifurcations = compute_bif_diag()
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :blue, :red), linewidth = 3.0)
        scatter!(ax1, bifurcations[1,1], bifurcations[1,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)

        # Export the figure
        savefig("ramped_SN_ews.pdf", fig)
end

# Execute the main
main()
