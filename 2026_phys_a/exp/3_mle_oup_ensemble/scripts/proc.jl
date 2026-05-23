"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""


# Import all the exported solutions
function analyse(case::Integer, θ)
        # Define the counting index
        start = 1

        # Loop over the number of executions of the main algorithm
        coefficients = Matrix{Float64}(undef, convert(Integer, Ne)*N_exec, 3)
        for index in 1:N_exec
                # Import solutions
                solutions = readin("α=$(α[case])/$index.csv")

                # Define the bounding index
                finish = (start - 1) + convert(Integer, Ne)

                # Loop over the columns of the imported solution matrix
                for (col_idx, col) in enumerate(eachcol(solutions))
                        # Compute the bounding index 

                        # Append the column in the larger matrix
                        coefficients[start:finish, col_idx] = col
                end

                # Update the counting index
                start = finish + 1
        end

        # Loop over the estimated coefficients
        for (index, column) in enumerate(eachcol(coefficients))
                # Compute the (linear) plotting index
                plt_idx = (case - 1)*3 + index 

                #=
                println("
----------------------- Column $(index): \n
                        Minimum                 (0%) = $(minimum(column)) 
                        Fifth percentile        (5%) = $(quantile(column, 0.05))
                        First quartile         (25%) = $(quantile(column, 0.25))
                        Median                 (50%) = $(median(column))
                        Third quartile         (75%) = $(quantile(column, 0.75))
                        Ninetyfifth percentile (95%) = $(quantile(column, 0.95))
                        Maximum               (100%) = $(maximum(column))
                       ")
                =#

                # Remove outliers from the solutions (consider the 1%~99% range)
                central = filter(x -> quantile(column, 0.01) ≤ x ≤ quantile(column, 0.99), column)

                # Build and plot the histogram of the solutions
                bins, pdf = fit_distribution(central, n_bins=Nb)
                barplot!(axes[plt_idx], bins, pdf, color = (:red,0.5), strokecolor = :black, strokewidth = 2)
                lines!(axes[plt_idx], [θ[index], θ[index]], [0, 1.1*maximum(pdf)], color = :blue, linewidth = 5.0)

                # Format the axis
                (axes[plt_idx]).title = L"\text{standard error} = %$(trunc(std(central)/sqrt(length(central)); digits=6))"
                (axes[plt_idx]).limits = (nothing, nothing, 0, 1.1*maximum(pdf))
                (axes[plt_idx]).xticks = [minimum(bins), θ[index], maximum(bins)] 
                (axes[plt_idx]).xtickformat = values -> ["$(trunc(value, digits=3))" for value in values]
                (axes[plt_idx]).yticks = [0, maximum(pdf)] 
                (axes[plt_idx]).ytickformat = values -> ["$(trunc(value, digits=0))" for value in values]
        end
end
