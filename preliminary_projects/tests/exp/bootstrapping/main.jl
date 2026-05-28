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
        # Loop over the parameter values
        for (parameter_index, μ) in enumerate(μ_set[1])
                # Solve the ensemble problem 
                x0 = [sqrt(-μ), μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the regularization coefficient values
                for (reg_index, α) in enumerate(α_set[1])
                        # Generate the figure
                        include("./scripts/figs.jl")

                        # Define the coefficients of the ground truth potential and plot it
                        β = [μ, 0.0, 1.0/3.0]
                        domain = LinRange(-2,2,1000)
                        lines!(ax1, domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                        ax1.ylabel = L"\text{V(x;\mu = %$(μ))}"

                        # Loop over the ensemble solutions
                        K = Vector{Float64}(undef, convert(Integer, Ne))
                        θ = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                        printstyled("μ = $(μ): solving the LLS problems\n"; bold=true, underline=true, color=:light_blue)
                        @showprogress for (solution_index, solution) in enumerate(ensemble.state)
                                # Check for tipping
                                if length(solution) == length(ensemble.time) 
                                        # Solve the LLS problem
                                        output = solve_lls(solution, α=α)
                                        θ[solution_index,:] = output.solution 
                                        K[solution_index] = output.cond_number

                                        # Plot and export the local information
                                        Fig = Figure()
                                        Ax = Axis(Fig[1,1], title = L"cond. number = %$(trunc(K[solution_index], digits=3))")
                                        lines!(Ax, domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                                        lines!(Ax, domain, [V(x, θ[solution_index,:]) for x in domain], color = :red, linewidth = 3.0)
                                        savefig("solutions/$solution_index.png", Fig)
                                else
                                        # Interrupt the execution and throw an error
                                        throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                                end
                        end

                        # Plot the ensemble mean potential and its std
                        θ_mean = vec(mean(θ, dims=1))
                        θ_std = vec(std(θ, dims=1))
                        display(θ_std)
                        lines!(ax1, domain, [V(x, θ_mean) for x in domain], color = (:red,0.25), linewidth = 3.0)

                        # Plot the histogram of the condition number
                        bins, pdf = fit_distribution(K, n_bins=Nb)
                        barplot!(ax2, bins, pdf, color = (:red,0.5), strokecolor = :black, strokewidth = 2)

                        # Export the figure
                        #savefig("bootostrapping/$(parameter_index)/$(reg_index).png", fig)
                        savefig("bootostrapping.pdf", fig)
                end
        end
end

# Execute the main
main()
