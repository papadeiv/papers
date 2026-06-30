"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/figs.jl")
include("./scripts/proc.jl")

# Define the main algorithm
function main()
        # Loop over the parameter values
        for (parameter_index, μ) in enumerate(μ_set[end])
                # Plot the ground truth in the search space
                scatter!(ax1, Point3f([μ, 0.0, 1.0/3.0]), markersize = 50, color = :blue, strokewidth = 1.0)

                # Solve the ensemble problem 
                x0 = [sqrt(-μ), μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the regularization values
                printstyled("μ = $(μ): solving the regularized LLS problems\n"; bold=true, underline=true, color=:light_blue)
                θ_mean = Matrix{Float64}(undef, length(α_set), 3)
                @showprogress for (regularization_index, α) in enumerate(α_set)
                        # Loop over the ensemble solutions
                        θ = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                        @showprogress for (solution_index, solution) in enumerate(ensemble.state)
                                # Check for tipping
                                if length(solution) == length(ensemble.time) 
                                        # Compute the regularized LLS solutions 
                                        θ[solution_index,:] = solve_lls(solution, α=α)
                                else
                                        # Interrupt the execution and throw an error
                                        throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                                end
                        end

                        # Compute the ensemble mean solution
                        θ_mean[regularization_index,:] = vec(mean(θ, dims=1))

                        # Plot the RLLS solutions in the search space
                        scatter!(ax1, [Point3f(θ_star) for θ_star in eachrow(θ)], markersize = 20, 
                                 color = regularization_index, colormap = (:Accent_6,0.05), colorrange = (1,length(α_set)))
                        if regularization_index > 1 && regularization_index < 3
                                scatter!(ax2, [Point3f(θ_star) for θ_star in eachrow(θ)], markersize = 20, 
                                 color = regularization_index, colormap = (:Accent_6,0.25), colorrange = (1,length(α_set)))
                        elseif regularization_index ≥ 3
                                scatter!(ax2, [Point3f(θ_star) for θ_star in eachrow(θ)], markersize = 20, 
                                 color = regularization_index, colormap = (:Accent_6,0.25), colorrange = (1,length(α_set)))
                                scatter!(ax3, [Point3f(θ_star) for θ_star in eachrow(θ)], markersize = 20, 
                                 color = regularization_index, colormap = (:Accent_6,1), colorrange = (1,length(α_set)))
                        end
                end

                # Plot the trajectory of the ensemble mean solution in the search space at varying regularization parameter
                θ_mean = vcat(transpose([μ, 0, 1.0/3.0]), θ_mean)
                lines!(ax1, θ_mean[1:end-1,1], θ_mean[1:end-1,2], θ_mean[1:end-1,3], color = :black, linewidth = 2.0)
                lines!(ax2, θ_mean[1:end-1,1], θ_mean[1:end-1,2], θ_mean[1:end-1,3], color = :black, linewidth = 2.0)
                for index in 1:length(α_set)
                        scatter!(ax1, θ_mean[index,1], θ_mean[index,2], θ_mean[index,3], markersize = 10, strokewidth = 1.0,
                                 color = index, colormap = :Accent_6, colorrange = (1,length(α_set)))
                        scatter!(ax2, θ_mean[index,1], θ_mean[index,2], θ_mean[index,3], markersize = 10, strokewidth = 1.0,
                                 color = index, colormap = :Accent_6, colorrange = (1,length(α_set)))
                end
                scatter!(ax2, Point3f([μ, 0.0, 1.0/3.0]), markersize = 30, color = :blue, strokewidth = 1.0)
        end

        # Export the figure
        savefig("tikhonov_regularization.pdf", fig)
end

# Execute the main
main()
