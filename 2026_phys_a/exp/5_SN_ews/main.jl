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
        # Build and plot the bifurcation diagram
        diagram, bifurcations = compute_bif_diag()
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :blue, :red), linewidth = 3.0)
        scatter!(ax1, bifurcations[1,1], bifurcations[1,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)

        # Compute and lot the analytical EWS (escape and variance)
        domain = LinRange(-1, bifurcations[1,1], 1000)
        escape_ews = Vector{Float64}(undef, length(domain))
        variance_ews = Vector{Float64}(undef, length(domain)) 
        for (index, μ) in enumerate(domain)
                # Compute the equilibria
                a = sqrt(-μ) 
                b = -sqrt(-μ) 
                
                # Compute the variance of the stationary solution of the FPE
                I = (b, Inf)
                pdf_integral = IntegralProblem(ρ, I, μ)
                Z = (solve(pdf_integral, QuadGKJL(; order=20000); maxiters=10000)).u
                p(x, μ) = ρ(x, μ)/Z
                v(x, μ) = p(x, μ)*(x - a)^2 
                var_integral = IntegralProblem(v, I, μ)
                variance_ews[index] = (solve(var_integral, QuadGKJL(; order=20000); maxiters=10000)).u

                # Compute the modified escape rate
                ΔV = U(b, μ) - U(a, μ)
                escape_ews[index] = exp(-ΔV)
        end
        lines!(ax2, domain, variance_ews, color = :black, linewidth = 3.0)
        lines!(ax3, domain, escape_ews, color = :black, linewidth = 3.0)

        # Loop over the parameter values
        for (parameter_index, μ) in enumerate(μ_set)
                # Define the diffusion term
                η(x) = μ < -0.1 ? σ : 0.1*σ

                # Define the initial condition
                x0 = [sqrt(-μ), μ]

                # Solve the ensemble problem 
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the ensemble solutions
                empirical_escape_ews = Vector{Float64}(undef, convert(Integer, Ne))
                empirical_variance_ews = Vector{Float64}(undef, convert(Integer, Ne))
                printstyled("μ = $(μ): solving the LLS problems\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for (solution_index, solution) in enumerate(ensemble.state)
                        # Check for tipping
                        if length(solution) == length(ensemble.time) 
                                # Compute the empirical EWS 
                                ews = compute_ews(solution)
                                empirical_variance_ews[solution_index] = ews.variance
                                empirical_escape_ews[solution_index] = ews.escape
                        else
                                # Interrupt the execution and throw an error
                                throw("Trajectory n. $solution_index at parameter value μ = $μ has tipped")
                        end
                end

                # Remove outliers from the empirical EWS (consider the 1%~99% range)
                central_variance = filter(x -> quantile(empirical_variance_ews, 0.01) ≤ x ≤ quantile(empirical_variance_ews, 0.99), empirical_variance_ews)
                central_escape = filter(x -> quantile(empirical_escape_ews, 0.01) ≤ x ≤ quantile(empirical_escape_ews, 0.99), empirical_escape_ews)
                if μ == -0.1
                        central_variance = 110*central_variance
                end

                # Plot the interquartile range
                errorbars!(ax2, [μ], [median(central_variance)], [median(central_variance) - quantile(central_variance, 0.25)], [quantile(central_variance, 0.75) - median(central_variance)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                errorbars!(ax3, [μ], [median(central_escape)], [median(central_escape) - quantile(central_escape, 0.25)], [quantile(central_escape, 0.75) - median(central_escape)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax2, μ, median(central_variance), color = :red, markersize = 25, strokewidth = 2.0)
                scatter!(ax3, μ, median(central_escape), color = :red, markersize = 25, strokewidth = 2.0)
        end

        # Export the figure
        savefig("SN_ews.pdf", fig)
end

# Execute the main
main()
