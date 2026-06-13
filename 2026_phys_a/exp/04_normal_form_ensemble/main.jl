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
        for (parameter_index, μ) in enumerate(μ_set)
                # Define the coefficients of the ground truth potential and plot it
                β = [μ, 0.0, 1.0/3.0]
                domain = LinRange(-2,2,1000)
                lines!((axes[parameter_index])[4], domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                ((axes[parameter_index])[4]).ylabel = L"\text{V(x;\mu = %$(μ))}"

                # Define the initial condition
                x0 = [sqrt(-μ), μ]

                # Solve the ensemble problem 
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the ensemble solutions
                θ = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                printstyled("μ = $(μ): solving the LLS problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for (solution_index, solution) in enumerate(ensemble.state)
                        # Check for tipping
                        if length(solution) == length(ensemble.time) 
                                # Solve the LLS problem
                                θ[solution_index,:] = solve_lls(solution)
                        else
                                # Interrupt the execution and throw an error
                                throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                        end
                end

                # Loop over the columns of the solution matrix
                mean_θ = Vector{Float64}(undef, 3)
                for (column_index, column) in enumerate(eachcol(θ))
                        # Remove outliers from the solutions (consider the 1%~99% range)
                        central = filter(x -> quantile(column, 0.01) ≤ x ≤ quantile(column, 0.99), column)

                        # Compute the mean of the solutions
                        mean_θ[column_index] = mean(central)

                        # Build and plot the histogram of the solutions
                        bins, pdf = fit_distribution(central, n_bins=Nb)
                        barplot!((axes[parameter_index])[column_index], bins, pdf, color = (:red,0.5), strokecolor = :black, strokewidth = 2)
                        lines!((axes[parameter_index])[column_index], [β[column_index], β[column_index]], [0, 1.1*maximum(pdf)], color = :blue, linewidth = 5.0)

                        # Format the axis
                        ((axes[parameter_index])[column_index]).limits = (nothing, nothing, 0, 1.1*maximum(pdf))
                        ((axes[parameter_index])[column_index]).xticks = [minimum(bins), β[column_index], maximum(bins)] 
                        ((axes[parameter_index])[column_index]).xtickformat = values -> ["$(trunc(value, digits=1))" for value in values]
                        ((axes[parameter_index])[column_index]).yticks = [0, maximum(pdf)] 
                        ((axes[parameter_index])[column_index]).ytickformat = values -> ["$(trunc(value, digits=1))" for value in values]
                end

                # Plot the ensemble mean approximation potential
                lines!((axes[parameter_index])[4], domain, [V(x, mean_θ) for x in domain], color = (:red, 0.5), linewidth = 3.0)
        end

        # Export the figure
        savefig("lls_solutions.pdf", fig)
end

# Execute the main
main()
