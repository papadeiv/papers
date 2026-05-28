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
        for (parameter_index, μ) in enumerate(μ_set)
                # Solve the ensemble problem 
                x0 = [sqrt(-μ), μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the regularization coefficient values
                printstyled("μ = $(μ): solving the RLLS problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for (reg_index, α) in enumerate(α_set)
                        # Generate the figure
                        include("./scripts/figs.jl")

                        # Define the coefficients of the ground truth potential and plot it
                        β = [μ, 0.0, 1.0/3.0]
                        domain = LinRange(-2,2,1000)
                        lines!(axes[6], domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                        (axes[6]).ylabel = L"\text{V(x;\mu = %$(μ))}"

                        # Loop over the ensemble solutions
                        θ = Matrix{Float64}(undef, convert(Integer, Ne), 3)
                        empirical_ews = Matrix{Float64}(undef, convert(Integer, Ne), 2)
                        for (solution_index, solution) in enumerate(ensemble.state)
                                # Check for tipping
                                if length(solution) == length(ensemble.time) 
                                        # Solve the LLS problem
                                        θ[solution_index,:] = solve_lls(solution, α=α)
                                        # Compute the empirical EWS 
                                        ews = compute_ews(solution, α=α)
                                        empirical_ews[solution_index,1] = ews.variance
                                        empirical_ews[solution_index,2] = ews.escape
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
                                barplot!(axes[column_index], bins, pdf, color = reg_index, colormap = (:phase,0.5), colorrange = (1,length(α_set)), strokecolor = :black, strokewidth = 2)
                                lines!(axes[column_index], [β[column_index], β[column_index]], [0, 1.1*maximum(pdf)], color = :blue, linewidth = 5.0)

                                # Format the axis
                                (axes[column_index]).limits = (nothing, nothing, 0, 1.1*maximum(pdf))
                                (axes[column_index]).xticks = [minimum(bins), β[column_index], maximum(bins)] 
                                (axes[column_index]).xtickformat = values -> ["$(trunc(value, digits=1))" for value in values]
                                (axes[column_index]).yticks = [0, maximum(pdf)] 
                                (axes[column_index]).ytickformat = values -> ["$(trunc(value, digits=1))" for value in values]
                        end

                        # Plot the ensemble mean approximation potential
                        lines!(axes[6], domain, [V(x, mean_θ) for x in domain], color = reg_index, colormap = (:phase,0.5), colorrange = (1,length(α_set)), linewidth = 3.0)

                        # Compute the analytcal EWS
                        a = sqrt(-μ) 
                        b = -sqrt(-μ) 
                        I = (b, Inf)
                        pdf_integral = IntegralProblem(ρ, I, μ)
                        Z = (solve(pdf_integral, QuadGKJL(; order=20000); maxiters=10000)).u
                        p(x, μ) = ρ(x, μ)/Z
                        v(x, μ) = p(x, μ)*(x - a)^2 
                        var_integral = IntegralProblem(v, I, μ)
                        true_variance = (solve(var_integral, QuadGKJL(; order=20000); maxiters=10000)).u
                        ΔV = U(b, μ) - U(a, μ)
                        true_escape = exp(-ΔV)
                        analytical_ews = [true_variance, true_escape]

                        # Loop over the columns of the EWS matrix
                        for (column_index, column) in enumerate(eachcol(empirical_ews))
                                # Remove outliers from the solutions (consider the 1%~99% range)
                                central = filter(x -> quantile(column, 0.01) ≤ x ≤ quantile(column, 0.99), column)

                                # Build and plot the histogram of the solutions
                                bins, pdf = fit_distribution(central, n_bins=Nb)
                                barplot!(axes[column_index+3], bins, pdf, color = reg_index, colormap = (:phase,0.5), colorrange = (1,length(α_set)), strokecolor = :black, strokewidth = 2)
                                lines!(axes[column_index+3], [analytical_ews[column_index], analytical_ews[column_index]], [0, 1.1*maximum(pdf)], color = :blue, linewidth = 5.0)

                                # Format the axis
                                (axes[column_index+3]).limits = (nothing, nothing, 0, 1.1*maximum(pdf))
                                (axes[column_index+3]).xticks = [minimum(bins), analytical_ews[column_index], maximum(bins)] 
                                (axes[column_index+3]).xtickformat = values -> ["$(trunc(value, digits=4))" for value in values]
                                (axes[column_index+3]).yticks = [0, maximum(pdf)] 
                                (axes[column_index+3]).ytickformat = values -> ["$(trunc(value, digits=0))" for value in values]
                        end

                        # Plot the current value of the regularisation coefficient on the colorbar and export the figure 
                        lines!(colorbar, [α, α], [0, 1], color = :black, linewidth = 3.0)
                        savefig("regularization_stationary/$(parameter_index)/$(reg_index).png", fig)
                end
        end
end

# Execute the main
main()
