function Stommel()
        # System parameters
        η1 = 3.00                                         # Temperature meridional difference 
        η3 = 0.30                                         # Timescale ratio between temperature and salinity 
        μ0 = 0.00                                         # Initial value of the bifurcation parameter
        μf = 1.20                                         # Final value of the bifurcation parameter
        ε = 5e-5                                          # Timescale separation
        σ = 0.050                                         # Noise level (additive)
        D = (σ^2)/2.0                                     # Diffusion level (additive) 

        # Dynamical system  
        f1(x, y, μ) = η1 - x - abs(x - y)*x               # Temperature (drift)
        f2(x, y, μ) = μ - η3*y - abs(x - y)*y             # Salinity (drift)
        f = [f1, f2]                                      # Drift
        Λ(t) = ε                                          # Shift/Ramp
        g(x, y) = σ                                       # Noise level (additive) 
        η = [g, g]                                        # Diffusion

        # Simulation parameters
        dt = 1e-1                                         # Timestep
        Ne = 2e2                                          # Number of particles in the ensemble 

        # Solve the ensemble problem 
        z0 = [get_equilibria(f1, f2, μ0).stable[1]; μ0]
        ensemble = evolve(f, η, Λ, z0, endparameter=μf, stepsize=dt, particles=Ne)
        t = ensemble.time
        μ = ensemble.parameter

        # Loop over the ensemble's sample paths
        μ_min, μ_max = μ0, μf 
        for (solution_index, solution) in enumerate(ensemble.state)
                # Extract the solution's coordinates
                x = solution[:,1]
                y = solution[:,2]

                # Construct the scalar observable
                ψ = x - y 

                # Identify the tipping and truncate the solution up to the that point
                tipping_index = find_tipping(ψ, 0.5)
                μt = μ[1:tipping_index]
                ψt = ψ[1:tipping_index]

                # Compute and plot the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(ψt; alg = "emd", n_modes=1).residuals

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

                # Loop over the window strides 
                ews = Matrix{Float64}(undef, length(subseries.trajectories), 5)
                printstyled("Stommel test case: progressing trajectory $solution_index of $(convert(Integer, Ne))\r"; bold=true, underline=true, color=:light_magenta)
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
                writeout(ews, "ews/stommel/$solution_index.csv")
        end

        # Compute the number of steps in the filtered subseries
        ews = readin("ews/stommel/1.csv")
        mask = (ews[:, 1] .>= μ_min) .& (ews[:, 1] .<= μ_max)
        filtered_ews = ews[mask, :]
        Nt_filtered = size(filtered_ews, 1)
        μ_filtered = filtered_ews[:,1]

        # Loop over the trajectories
        ensemble_θ1 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_θ2 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_θ3 = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        ensemble_variance = Matrix{Float64}(undef, convert(Integer, Ne), Nt_filtered)
        printstyled("\n  Computing ensemble means\n"; bold=true, underline=true, color=:light_green)
        @showprogress for n in 1:convert(Integer, Ne)
                # Import the matrix from file
                ews = readin("ews/stommel/$n.csv")

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
        lines!(ax1R, μ_filtered, ensemble_mean_variance, color = (:red,0.75), linewidth = 3.0)
        lines!(ax2R, μ_filtered, ensemble_mean_ews, color = (:red,0.75), linewidth = 3.0)

        # Plot the tipping point and format the axes
        lines!(ax1R, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        lines!(ax2R, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        println(ensemble_mean_variance[end])
        println(ensemble_mean_ews[end])

        # Export the figure
        savefig("ews_comparison.pdf", fig)
end
