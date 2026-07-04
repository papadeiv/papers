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
        # Find the minimum interval
        μ_min, μ_max = 0.0, 1.2
        for n in 1:convert(Integer, Ne)
                # Import the matrix from file
                ews = readin("ews/stommel/$n.csv")
                #println(n, ") μ_min=", ews[1,1], " (", μ_min, "), ", ews[end,1], " (", μ_max, ")")

                if ews[1,1] > μ_min
                        #println("------> Update!")
                        μ_min = ews[1,1]
                end
                if ews[end,1] < μ_max
                        #println("------> Update!")
                        μ_max = ews[end,1]
                end
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
        printstyled("\nComputing ensemble means\n"; bold=true, underline=true, color=:light_magenta)
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
        lines!(ax1, μ_filtered, ensemble_mean_variance, color = (:red,0.75), linewidth = 3.0)
        lines!(ax2, μ_filtered, ensemble_mean_ews, color = (:red,0.75), linewidth = 3.0)

        # Plot the tipping point and format the axes
        lines!(ax1, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        lines!(ax2, [μ_max, μ_max], [-2,2], color = :black, linewidth = 3.0, linestyle = :dash)
        ax1.limits = (μ_min, μf, 0, 0.003)
        ax1.xticks = [μ_min, μf]
        ax1.yticks = [0, 0.003]
        ax2.limits = (μ_min, μf, 0, 1)
        ax2.xticks = [μ_min, μf]
        ax2.yticks = [0,1]

        # Export the figure
        savefig("ramped_stommel_ews.pdf", fig)
end

# Execute the main
main()
