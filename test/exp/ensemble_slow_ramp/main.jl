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
        # Solve the ensemble problem 
        sample_paths = evolve_ramped_1d(f, g, η, x0, δt=δt, μf=μf, Ne=Ne)
        t = sample_paths.time
        μ = sample_paths.parameter

        # Plot the true scalar potential
        plot_ground_truth(μ)

        # Compute the quasi-steady equilibrium in the windowed parameter range of the subseries 
        qse = [(get_equilibria(f, μc, domain=[-10,10])).stable[2] for μc in μ] 

        # Loop over the ensemble's sample paths
        printstyled("Computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:Ne
                # Extract the current solution from the ensemble and detrend it
                u = sample_paths.states[n]
                u_det = detrend(u, qse = qse)
                u_res = u_det.residuals

                # Solve the NLLS problem
                c[n,:] = fit_potential(u_res, n_coeff=Nc, noise=σ, optimiser=β, attempts=Na)

                # Reconstruct a shifted potential (to match the ground truth)
                xs, Vs = shift_potential(U, x0[1], μ[end], c[n,:])

                # Compute the random variables targeted by the analysis 
                parameters[n,1] = xs                    # Estimated stable equilibrium
                parameters[n,2] = V(xs,c[n,:])          # Potential value at equilibrium
                parameters[n,3] = Vxx(xs,c[n,:])        # Curvature value at equilibrium
                parameters[n,4] = exp(V(xs,c[n,:])/D)   # Escape probability at equilibrium

                # Plot the results (only for the first 100 particles in the ensemble) 
                if n < 100
                        plot_solutions(μ, u, u_res)       # Timeseries of the ensemble
                        plot_reconstruction(Vs)       # Timeseries of the ensemble
                end
        end

        # Plot the drift of the quasi-steady equilibrium (QSE) and print the info
        plot_drift(μ, qse)
        print_info(length(t))

        # Plot the random variables
        plot_coeffs(c)                                  # Solutions of the NLLS problem
        plot_rvs(parameters)                            # Transformations of the coefficients 

        # Export the figures
        savefig("ensemble_slow_ramp/solutions.png", fig1)
        savefig("ensemble_slow_ramp/analysis.png", fig4)
end

# Execute the main
main()
