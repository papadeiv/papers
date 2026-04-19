"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")

# Define the main algorithm
function main()
        #--------------------#
        #        Wood        #
        #--------------------#
        
        # Solve the ensemble slow-fast SDEs
        equilibria = get_equilibria(f1, f2, H0, guesses=[[Sn0,St0], [-0.5,0.5], [-0.1,0.6]])
        z0 = [equilibria.stable[1]; H0]
        ensemble = evolve(f, η, Λ, z0, stepsize=δt, steps=Nt)

        # Extract the x and y components
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]

        # Extract the parameter's shift and timesteps
        H = ensemble.parameter
        t = ensemble.time

        # Plot the timeseries and phase space trajectory
        fig1 = Figure()
        fig2 = Figure()
        ax1 = Axis(fig1[1,1], xlabel = L"S_N", ylabel = L"S_T", title = "Wood's 5-box model'")
        lines!(ax1, x, y, color = (:magenta,0.35), linewidth = 0.5)
        ax2 = Axis(fig2[1,1], ylabel = L"S_N", title = "Wood's 5-box model'")
        lines!(ax2, t, x, color = :red, linewidth = 1.0)
        ax3 = Axis(fig2[2,1], xlabel = L"t", ylabel = L"S_T")
        lines!(ax3, t, y, color = :blue, linewidth = 1.0)

        #-------------------#
        #    Saddle-node    #
        #-------------------#
        
        # Solve the ensemble slow-fast SDEs
        equilibria = get_equilibria(F1, F2, μ0)
        z0 = [equilibria.stable[1]; μ0]
        ensemble = evolve(F, ζ, Λ, z0, stepsize=δt, steps=Nt)

        # Extract the x and y components
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]

        # Extract the parameter's shift and timesteps
        M = ensemble.parameter
        t = ensemble.time

        ax4 = Axis(fig1[1,2], xlabel = L"x", ylabel = L"y", title = "Saddle-node normal form")
        lines!(ax4, x, y, color = (:magenta,0.35), linewidth = 0.5)
        ax5 = Axis(fig2[1,2], ylabel = L"x", title = "Saddle-node normal form")
        lines!(ax5, t, x, color = :red, linewidth = 1.0)
        ax6 = Axis(fig2[2,2], xlabel = L"t", ylabel = L"y")
        lines!(ax6, t, y, color = :blue, linewidth = 1.0)

        savefig("comparison/phase_space.png", fig1)
        savefig("comparison/timeseries.png", fig2)

        #display(Υ/Vt)
        #display(H[end])
        #display(M[end])
end

# Execute the main
main()
