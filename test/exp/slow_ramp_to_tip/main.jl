"""
    Escape probability EWS to a slow-fast saddle-node.

We test the performance of the escape probability EWS computed out of the reconstructed potential in action.
By that we mean that we apply it to a synthetic timeseries that replicates the intended target of the method:
- it is the solution of a SDE with non-autonomous drift term, meaning that the bifurcation parameter is not stationary but instead it has a time dependence;
- the parameter shift or ramp is linear and with a very small slope w.r.t. time meaning that the system undergoes slow passage thorugh a saddle-nome bifurcation generating the tipping.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Solve the slow-fast SDE 
        sample_path = evolve_ramped_1d(f, g, η, x0, δt=δt, μf=μf)
        t = sample_path.time
        μ = sample_path.parameter
        u = (sample_path.states)[1]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(μ, u, width)
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Loop over the ensemble's sample paths
        printstyled("Computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:Ne
                # Compute the quasi-steady equilibrium in the windowed parameter range of the subseries 
                qse = [(get_equilibria(f, μc, domain=[-10,10])).stable[2] for μc in ensemble.timesteps[n]] 
                
                # Extract the current solution from the ensemble and detrend it 
                detrended_solution = detrend(ensemble.trajectories[n], qse = qse)
                local u = detrended_solution.residuals

                # Solve the NLLS problem
                coefficients = fit_potential(u, n_coeff=Nc, noise=σ, optimiser=β, attempts=Na)
                push!(solutions, coefficients)

                # Reconstruct a shifted potential (to match the ground truth)
                xs, Vs = shift_potential(U, x0[1], ensemble.timesteps[n][end], coefficients)

                # Compute the random variables targeted by the analysis 
                push!(parameters, [xs,                          # Estimated stable equilibrium
                                   V(xs,coefficients),          # Potential value at equilibrium
                                   Vxx(xs,coefficients),        # Curvature value at equilibrium
                                   exp(V(xs,coefficients)/D)])  # Escape probability at equilibrium

                # Get local extrema of the true and reconstructed potential and redefine their functions
                bounds = get_bounds(ensemble.timesteps[n][end], coefficients)
                Ur(x) = U(x, ensemble.timesteps[n][end])
                Urxx(x) = Uxx(x, ensemble.timesteps[n][end])
                Vr(x) = V(x, coefficients)
                Vrxx(x) = Vxx(x, coefficients)

                # Compute the escape EWS
                ground_truth = estimate_escape(Ur, Urxx, bounds.true_min, bounds.true_max, σ)
                reconstruction = estimate_escape(Vr, Vrxx, bounds.approx_min, bounds.approx_max, σ)
                push!(escapes, [ground_truth.prefactor,         # Ground truth prefactor 
                                ground_truth.LDP,               # Ground truth large deviation 
                                reconstruction.prefactor,       # Reconstruction prefactor
                                reconstruction.LDP,             # Reconstruction large deviation 
                               ])

                # Compute the reconstruction error
                push!(error, get_error(ensemble.timesteps[n][end], Vs))

                # Plot and export the residuals, distributions and potentials
                if mod(n,100) == 0
                        fig2 = plot_solutions(ensemble.timesteps[n], 
                                      ensemble.trajectories[n], 
                                      detrended_solution, 
                                      coefficients, 
                                      Vs)
                        savefig("slow_ramp_to_tip/solutions/$n.png", fig2)
                end
                
        end

        # Plot the results
        plot_tip(t, u, ensemble.tipping_point)                  # Non-stationary sample path
        plot_coeffs(ensemble.timesteps, solutions, error)       # Timeseries of the solutions and numerical error 
        plot_rvs(parameters)                                    # Transformations of the coefficients
        plot_ews(ensemble.timesteps, escapes)                   # Escape probability EWS 

        # Export the figures
        savefig("slow_ramp_to_tip/timeseries.png", fig1)
        savefig("slow_ramp_to_tip/coefficients.png", fig3)
        savefig("slow_ramp_to_tip/analysis.png", fig7)
        savefig("slow_ramp_to_tip/components.png", fig11)
        savefig("slow_ramp_to_tip/ews.png", fig13)
end

# Execute the main
main()
