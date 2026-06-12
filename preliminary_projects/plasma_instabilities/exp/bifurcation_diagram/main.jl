# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/figs.jl")

# Define the main algorithm
function main()
        # Compute and plot the bifurcation diagram
        equilibria = get_equilibria(f, μ_set[1], domain=[-10,10])
        diagram, bifurcations = compute_bif_diag(equilibria.stable)
        lines!(ax, diagram[:,1], diagram[:,2], color = ifelse.(diagram[:,3] .≤ 0, :black, :gray), linewidth = 3.0)
 
        # Iterate over the parameters grid 
        for (parameter_index, μ) in enumerate(μ_set)
                # Solve the differential equation 
                solution = evolve(f, η, Λ, [x0, μ], stepsize=dt, steps=Nt)
                p = solution.parameter
                x = solution.state[1]

                # Plot the solution in the bifurcation diagram
                lines!(ax, p, x, color = parameter_index, colormap = (:phase,1.0), colorrange = (1,length(μ_set)), linewidth = 3.0)
                scatter!(ax, p[1], x[1], color = :blue, markersize = 20, strokewidth = 1.0)
                scatter!(ax, p[end], x[end], color = :red, markersize = 20, strokewidth = 1.0)
        end

        # Export the figure
        savefig("bifurcation_diagram.pdf", fig)
end

# Execute the main
main()
