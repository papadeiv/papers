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
        μ1 = bifurcations[2,1]
        μ2 = bifurcations[1,1]

        # Loop over the diffusion value
        subdomain = LinRange(μ1, μ2, 1000)
        for (diffusion_index, D) in enumerate(reverse(D_set))
                escape_ews = Vector{Float64}(undef, length(subdomain))
                prefactor = Vector{Float64}(undef, length(subdomain))
                kramers = Vector{Float64}(undef, length(subdomain))
                for (index, μ) in enumerate(subdomain)
                        # Compute the local minimum and maximum of the potential
                        equilibria = get_equilibria(f, μ, domain=[-10,10])
                        a = maximum(equilibria.stable)
                        b = maximum(equilibria.unstable)

                        # Compute the escape ews
                        ΔV = U(b, μ) - U(a, μ)
                        escape_ews[index] = exp(-ΔV/D)

                        # Compute Kramer's escape rate
                        prefactor[index] = sqrt(abs(Uxx(b, μ))*Uxx(a, μ))/(2*pi)
                        kramers[index] = prefactor[index]*escape_ews[index]
                end
                lines!(ax1, subdomain, escape_ews, color = diffusion_index, colormap = (:roma,1.0), colorrange = (1,length(D_set)), linewidth = 3.0)
                lines!(ax2, subdomain, prefactor, color = :black, linewidth = 3.0)
                lines!(ax3, subdomain, kramers, color = diffusion_index, colormap = (:roma,1.0), colorrange = (1,length(D_set)), linewidth = 3.0)
        end

        # Format and export the figure
        ax1.limits = (μ1, μ2, 0, 1)
        ax1.yticks = [0, 1]
        ax2.limits = (μ1, μ2, 0, 0.06)
        ax2.yticks = [0, 0.06]
        ax3.limits = (μ1, μ2, 0, 0.03)
        ax3.yticks = [0, 0.03]
        savefig("kramers_formula.pdf", fig)
end

# Execute the main
main()
