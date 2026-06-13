"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Compute the initial conditions for the continuation problem
        equilibria = get_equilibria(f1, f2, μ.η2)

        # Compute and plot the upper branch of the bifurcation diagram
        diagram, bifurcations[1,:] = compute_bif_diag((equilibria.stable)[2], [0.0,2.0])
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,4] .≤ 0, :blue, :red), linewidth = 3.0)
        lines!(ax2, diagram[:,1], diagram[:,3], color = ifelse.(diagram[:,4] .≤ 0, :blue, :red), linewidth = 3.0)

        # Compute and plot the lower branch of the bifurcation diagram
        diagram, bifurcations[2,:] = compute_bif_diag((equilibria.stable)[1], [0.9,2.0])
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,4] .≤ 0, :blue, :red), linewidth = 3.0)
        lines!(ax2, diagram[:,1], diagram[:,3], color = ifelse.(diagram[:,4] .≤ 0, :blue, :red), linewidth = 3.0)
        bifurcations[2,1] = 0.9
        bifurcations[2,3] = 3.0

        # Plot bounding box for the region of bistability
        poly!(ax1, Point2f[(bifurcations[2,1],-1), (bifurcations[1,1],-1), (bifurcations[1,1],1.5), (bifurcations[2,1],1.5)], color = (:black,0.05), strokecolor = :black, strokewidth = 1.0) 
        poly!(ax2, Point2f[(bifurcations[2,1],-1), (bifurcations[1,1],-1), (bifurcations[1,1],3.5), (bifurcations[2,1],3.5)], color = (:black,0.05), strokecolor = :black, strokewidth = 1.0) 

        # Plot the S-N bifurcations
        scatter!(ax1, bifurcations[1,1], bifurcations[1,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)
        scatter!(ax1, bifurcations[2,1], bifurcations[2,2], color = :green, markersize = 15, strokecolor = :black, strokewidth = 1.0)
        scatter!(ax2, bifurcations[1,1], bifurcations[1,3], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)
        scatter!(ax2, bifurcations[2,1], bifurcations[2,3], color = :green, markersize = 15, strokecolor = :black, strokewidth = 1.0)

        # Format the μ axes of both plots
        ax1.xticks = [0, bifurcations[1,1], bifurcations[2,1], 2]
        ax2.xticks = [0, bifurcations[1,1], bifurcations[2,1], 2]
        
        # Export the figure
        savefig("stommel_bif_diag.pdf", fig)
end

# Execute the main
main()
