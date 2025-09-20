"""
    ?

???
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Loop over the discrete ranges
        for m in 1:Nμ
                # Create empty layouts for the figures
                include("./scripts/figs.jl")

                # Solve the slowly ramped ensemble problem
                sample_paths = evolve_ramped_1d(f, g, η, x0[m,:], δt=δt, μf=μf[m], Ne=Ne)
                t = sample_paths.time
                μ = sample_paths.parameter

                # Plot the true scalar potential
                plot_ground_truth(μ)

                # Compute the quasi-steady equilibrium of the slow ramp 
                qse = [(get_equilibria(f, μc, domain=[-10,10])).stable[2] for μc in μ] 

                # Loop over the ensemble's sample paths
                printstyled("Range $m of $Nμ: ($(μ0[m]),$(μf[m]))\n"; bold=true, underline=true, color=:green)
                printstyled("Solving the least-squares problems across the ensemble\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:Ne
                        # Extract the current solution from the ensemble and detrend it
                        u = sample_paths.states[n]
                        u_det = detrend(u, qse = qse)
                        u_res = u_det.residuals

                        # Solve the NLLS problem
                        coefficients[n,:] = fit_potential(u_res, n_coeff=Nc, noise=σ, optimiser=β, attempts=Na)

                        # Compute the random variables for the analysis
                        push!(analysis, analyse(coefficients[n,:]))

                        # Reconstruct a shifted potential (to match the ground truth)
                        xs, Vs = shift_potential(U, x0[m,1], μ[end], coefficients[n,:])

                        # Compute the escape EWS
                        escape[n,m] = estimate_escape(coefficients[n,:], σ).LDP
              
                        # Plot the results (only for the first 100 particles in the ensemble) 
                        if n < 100
                                plot_solutions(μ, u, u_res)   # Timeseries of the ensemble
                                plot_reconstruction(Vs)       # Timeseries of the ensemble
                        end
                end

                # Plot the drift, info and results 
                plot_drift(μ, qse)
                print_info(length(t), m)
                plot_results(coefficients, analysis)

                # Export the figures
                savefig("ensemble_slow_ramp/solutions/$m.png", fig1)
                savefig("ensemble_slow_ramp/coefficients/$m.png", fig4)
                savefig("ensemble_slow_ramp/analysis/$m.png", fig7)
                savefig("ensemble_slow_ramp/ews/$m.png", fig15)
        end

        #=
        # Plot and export the escape ews
        true_escape = plot_ews(escape)
        savefig("ensemble_slow_ramp/ews.png", fig11)

        # Export data
        writeCSV(escape, "../../res/data/escape.csv")
        writeCSV(true_escape, "../../res/data/true_escape.csv")
        =#
end

# Execute the main
main()
