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
        # Build and plot the bifurcation diagram
        equilibria = get_equilibria(f, c0, domain=[-10,10])
        diagram, bifurcations = compute_bif_diag(equilibria.stable)
        lines!(ax1, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :blue, :red), linewidth = 3.0)

        # Loop over the slices of the parameter space
        domain = LinRange(0,9,1000)
        for (index, μ) in enumerate(c)
                # Plot the potential at the current slice
                lines!(axes[index], domain, [U(x,μ) for x in domain], color = :black, linewidth = 2.0)
        end

        # Compute and plot the escape ews
        subdomain = LinRange(bifurcations[2,1], bifurcations[1,1], 1000)
        escape_ews = Vector{Float64}(undef, length(subdomain))
        for (index, μ) in enumerate(subdomain)
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                a = maximum(equilibria.stable)
                b = maximum(equilibria.unstable)
                ΔV = U(b, μ) - U(a, μ)
                escape_ews[index] = exp(-ΔV)
        end
        lines!(ax3, subdomain, escape_ews, color = :black, linewidth = 3.0)

        # Plot bounding boxes for the region of bistability
        poly!(ax1, Point2f[(bifurcations[2,1],0), (bifurcations[1,1],0), (bifurcations[1,1],9), (bifurcations[2,1],9)], color = (:black,0.05), strokecolor = :black, strokewidth = 1.0) 
        poly!(ax3, Point2f[(bifurcations[2,1],0), (bifurcations[1,1],0), (bifurcations[1,1],1), (bifurcations[2,1],1)], color = (:black,0.05), strokecolor = :black, strokewidth = 1.0) 

        # Plot the S-N bifurcations
        scatter!(ax1, bifurcations[:,1], bifurcations[:,2], color = :yellow, markersize = 15, strokecolor = :black, strokewidth = 1.0)

        # Export the figure
        savefig("escape_ews.pdf", fig)
end

# Execute the main
main()
