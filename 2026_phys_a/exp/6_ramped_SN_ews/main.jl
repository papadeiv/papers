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
        #=
        # Compute and lot the analytical EWS (escape and variance)
        domain = LinRange(μ0, μf, 1000)
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
        lines!(ax3_L, domain, variance_ews, color = :black, linewidth = 3.0)
        lines!(ax3_R, domain, escape_ews, color = :black, linewidth = 3.0)
        =#

        # Solve the ensemble problem 
        x0 = [sqrt(-μ0), μ0]
        ensemble = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        t = ensemble.time
        μ = ensemble.parameter
        display(length(t))

        # Loop over the ensemble's sample paths
        for (solution_index, solution) in enumerate(ensemble.state)
                # Define the transparency index
                transparency = solution_index < convert(Integer, Ne) ? 0.05 : 1.0 

                # Extract the length of the sample path and truncate the parameter solution to match it
                N_end = length(solution)
                μ_end = μ[1:N_end]

                # Plot the sample path
                downsample_index = 40::Integer
                lines!(ax1, μ_end[1:downsample_index:end], solution[1:downsample_index:end], color = (:black,transparency), linewidth = 1.0)

                # Identify the tipping and truncate the solution up to the that point
                #tipping_index = find_tipping(solution; width = 300, threshold = 2.75, verbose=false) - 100::Integer
                tipping_index = find_tipping(solution; width = 100, threshold = 10.00, verbose=true) - 100::Integer
                μt = μ[1:tipping_index]
                xt = solution[1:tipping_index]

                # Compute and plot the quasi-stationary residuals of the truncated timeseries
                residuals = detrend(xt; alg = "emd", n_modes = 1).residuals
                lines!(ax2, μt[1:2*downsample_index:end], residuals[1:2*downsample_index:end], color = (:black, transparency), linewidth = 1.0)

                # Convert the sliding window across the residuals timeseries into an ensemble of suberies
                subseries = preprocess_solution(μt, residuals, window_size, verbose = false)

                # Loop over the window strides 
                ews = Matrix{Float64}(undef, length(subseries.trajectories), 3)
                printstyled("Progressing through the sliding window for trajectory $solution_index of $(convert(Integer, Ne))\n\r"; bold=true, underline=true, color=:light_magenta)
                @showprogress for (n, (parameter, trajectory)) in enumerate(zip(subseries.timesteps, subseries.trajectories))
                        #println("Trajectory $solution_index, parameter range = [$(parameter[1]), $(parameter[end])]")
                        # Compute the variance and escape early-warning signal
                        signal = compute_ews(trajectory, α=3e0)
                        ews[n,1] = parameter[end]
                        ews[n,2] = signal.variance 
                        ews[n,3] = signal.escape 
                        #println("-------------------------------------------")
                end

                # Plot the EWS 
                lines!(ax3_L, ews[:,1], ews[:,2], color = (:red,transparency), linewidth = 1.0)
                lines!(ax3_R, ews[:,1], ews[:,3], color = (:red,transparency), linewidth = 1.0)

                # Plot the tipping point of the last trajectory and format the axes
                if solution_index == convert(Integer, Ne)
                        lines!(ax2, [μt[end], μt[end]], [-1,1], color = :black, linewidth = 3.0, linestyle = :dash)
                        lines!(ax3_L, [μt[end], μt[end]], [-1,1], color = :black, linewidth = 3.0, linestyle = :dash)
                        lines!(ax3_R, [μt[end], μt[end]], [-1,1], color = :black, linewidth = 3.0, linestyle = :dash)
                        ax3_L.limits = (ews[1,1], 0, 0, 0.008)
                        ax3_L.xticks = [ews[1,1], 0]
                        ax3_L.yticks = [0, 0.008]
                        ax3_R.limits = (ews[1,1], 0, 0, 1)
                        ax3_R.xticks = [ews[1,1], 0]
                        ax3_R.yticks = [0, 1]
                end
        end

        # Build and plot the bifurcation diagram
        diagram, bifurcations = compute_bif_diag()
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :blue, :red), linewidth = 3.0)
        scatter!(ax1, bifurcations[1,1], bifurcations[1,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)

        # Export the figure
        savefig("stationary_SN_ews.pdf", fig)
end

# Execute the main
main()
