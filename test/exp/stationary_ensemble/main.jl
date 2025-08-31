# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Solve the ensemble problem 
        t, X = evolve_ensemble(f, η, μ, IC=[x0], δt=δt, Nt=Nt, Ne=Ne)

        # Loop over the ensemble's sample paths
        printstyled("Computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:Ne
                # Extract the current solution from the ensemble and center it
                u = X[n,:] .- x0 

                # Solve the NLLS problem
                c[n,:] = fit_potential(u, n_coeff=Nc, n_bins=Nb, noise=σ, optimiser=β, attempts=Na)

                # Reconstruct a shifted potential (to match the ground truth)
                xs, Vs = shift_potential(U, x0, μ, c[n,:])

                # Compute the random variables targeted by the analysis 
                parameters[n,1] = xs                    # Estimated stable equilibrium
                parameters[n,2] = V(xs,c[n,:])          # Potential value at equilibrium
                parameters[n,3] = Vxx(xs,c[n,:])        # Curvature value at equilibrium
                parameters[n,4] = exp(V(xs,c[n,:])/D)   # Escape probability at equilibrium

                # Compute the reconstruction error
                error[n] = get_error(μ, Vs)

                # Plot the results (only for the first 100 particles in the ensemble) 
                if n < 100
                        plot_ax1(t, u)                  # Timeseries of the ensemble
                        plot_ax2(Vs)                    # Reconstructed (shifted) potential 
                end
        end

        # Plot the random variables
        plot_coeffs(c)                                  # Solutions of the NLLS problem
        plot_error(error)                               # Numerical error of the reconstruction
        plot_rvs(parameters)                            # Transformations of the coefficients 

        # Export the data
        writeCSV(parameters, "../../res/data/parameters.csv")
        writeCSV(error, "../../res/data/error.csv")
        writeCSV(c, "../../res/data/coefficients.csv")

        # Export the figures
        savefig("solutions.png", fig1)
        savefig("analysis.png", fig3)
end

# Execute the main
main()
