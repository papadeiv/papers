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
        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/AMOC_transport_depth_0-1000m_monthly.nc")
        time = data["time"][:]
        observable = data["Transport"][:]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(time, observable, window_size)
        tipping = ensemble.tipping_point
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Loop over the ensemble's sample paths
        Threads.@threads for n in 1:100#Ne
                # Display the counter status
                c = atomic_add!(counter, 1)
                if c % 1 == 0
                        lock(print_lock) do
                                printstyled("\rComputing the least-squares solutions across the ensemble using $(Threads.nthreads()) threads: $c/$Ne"; bold=true, underline=true, color=:light_blue)
                                flush(stdout)
                        end
                end

                # Extract windowed subseries
                tw = copy(ensemble.timesteps[n])
                uw = copy(ensemble.trajectories[n]) 

                # Detrend it 
                detrended_solution = detrend(uw)
                residuals = detrended_solution.residuals

                # Solve the nonlinear least-squares problem to fit a cubic potential
                solution = fit_potential(residuals, noise=σ, transformation=[0.0,1.0,8.0], optimiser=β, attempts=Na)
                push!(solutions, solution.fit)

                # Perform postprocessing analysis on the solutions
                analyse(solution.fit)
                push!(time_idx, tw[end])
        end

        # Sort the ensemble solutions by temporal order
        sort_idx = sortperm(time_idx)
        time_idx = time_idx[sort_idx]
        solutions = solutions[sort_idx]
        ews = ews[sort_idx]

        # Loop over the solutions
        @showprogress for n in 1:100#Ne
                # Create empty layouts for the figures
                include("./scripts/figs.jl")

                # Plot the full timeseries
                lines!(ax1, time, observable, color = (:black,0.5), linewidth = 0.5)
                labels = ["0", "1", "2"]
                ax1.xtickformat = values -> ["$(label)" for label in labels]
                ax3.xtickformat = values -> ["$(label)" for label in labels]

                # Plot the sliding window
                poly!(ax1, Point2f[(time[n], -7), 
                                   (time_idx[n], -7),
                                   (time_idx[n], 25),
                                   (time[n], 25),
                                  ], color = (CtpGray, 0.1), strokecolor = :black, strokewidth = 1.0)

                # Plot the subseries
                lines!(ax1, [time[n], time_idx[n]], observable[n,Nt+n], color = (:black,1.0), linewidth = 0.5)

                # Plot the location of the tipping
                lines!(ax1, [time[tipping], time[tipping]], [-7, 25], color = (:red, 1.0), linewidth = 3.0, linestyle = :dash)
                lines!(ax3, [time[tipping], time[tipping]], [-0.5, 1.5], color = (:red, 1.0), linewidth = 3.0, linestyle = :dash)

                # Plot the ews
                lines!(ax3, time_idx[1:n], ews[1:n], color = CtpBlue, linewidth = 3.0)
                scatter!(ax3, time_idx[n], ews[n], color = CtpBlue, strokewidth = 3.0)

                # Plot the potential reconstruction
                domain = collect(LinRange(-5, 5, 1000))
                lines!(ax2, domain, [V(x) for x in domain], color = (:black, 1.00), linewidth = 3.0)

                # Export the figure
                savefig("amoc/$idx.png", fig1)
        end
end

# Execute the main
main()
