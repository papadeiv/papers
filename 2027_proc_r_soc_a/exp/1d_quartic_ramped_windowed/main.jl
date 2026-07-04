"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
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
        sample_path = evolve(f, η, Λ, x0, endparameter=μf, stepsize=δt)
        t = sample_path.time
        μ = sample_path.parameter
        u = (sample_path.state)[1]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(μ, u, width)
        tipping = ensemble.tipping_point
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        #=
        # Loop over the ensemble's sample paths
        Threads.@threads for n in 1:Ne
                # Display the counter status
                c = atomic_add!(counter, 1)
                if c % 1 == 0
                        lock(print_lock) do
                                printstyled("\rComputing the least-squares solutions across the ensemble using $(Threads.nthreads()) threads: $c/$Ne"; bold=true, underline=true, color=:light_blue)
                                flush(stdout)
                        end
                end

                # Extract windowed subseries
                μw = copy(ensemble.timesteps[n])
                uw = copy(ensemble.trajectories[n]) 

                # Compute the quasi-steady equilibrium in the windowed parameter range of the subseries 
                qse = [(get_equilibria(f, μc, domain=[-10,10])).stable[2] for μc in μw] 
                
                # Extract the current solution from the ensemble and detrend it 
                detrended_solution = detrend(uw, qse = qse)
                residuals = detrended_solution.residuals

                # Solve the nonlinear least-squares problem to fit a cubic potential
                solution = fit_potential(residuals, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                push!(solutions, [μw[end]; solution.fit])

                # Perform postprocessing analysis on the solutions
                Vs = analyse(solution.fit, μw[end])
        end

        # Export the postprocessed data and vacate the analysis arrays for the next iteration
        export_data(1)
        empty!(results)
        empty!(solutions)
        =#

        # Import the exported data
        data = import_data(1)

        #=
        # Loop over the solutions and results of the previous loop
        printstyled("\nExporting the movie for the potential reconstruction\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in union(1:100:Ne, [Ne])
                # Extract windowed subseries
                μw = copy(ensemble.timesteps[n])
                uw = copy(ensemble.trajectories[n]) 

                # Extract least-squares solutions
                coefficients = data.coefficients[n]

                # Compute shifted potential
                xs, Vs =  shift_potential(U, x0[1], μw[1], coefficients[2:4])

                # Compute the quasi-steady equilibrium in the windowed parameter range of the subseries 
                qse = [(get_equilibria(f, μc, domain=[-10,10])).stable[2] for μc in μw]
                
                # Extract the current solution from the ensemble and detrend it 
                detrended_solution = detrend(uw, qse = qse)
                residuals = detrended_solution.residuals

                # Plot and export the sample path, histogram and potential
                include("./scripts/figs.jl")
                plot_solution(μ, u, tipping, μw, residuals, Vs, n)
        end
        =#

        # Plot and export the ews timeseries
        plot_ews(μ, u, data.analysis, tipping)
end

# Execute the main
main()
