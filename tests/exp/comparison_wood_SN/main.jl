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
        # Create figures
        fig1 = Figure()
        fig2 = Figure()
        fig3 = Figure()
        fig4 = Figure()

        #--------------------#
        #        Wood        #
        #--------------------#
        
        # Solve the ensemble slow-fast SDE
        equilibria = get_equilibria(f1, f2, H0, guesses=[[Sn0,St0], [-0.5,0.5], [-0.1,0.6]])
        z0 = [equilibria.stable[1]; H0]
        ensemble = evolve(f, η, Λ, z0, stepsize=δt, steps=Nt)

        # Extract the solutions
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]
        H = ensemble.parameter
        t = ensemble.time
        T = range(0, 10, length(t))

        # Plot the timeseries and phase space trajectory
        ax1 = Axis(fig1[1,1], xlabel = L"S_N", ylabel = L"S_T", title = "Wood's 5-box model'")
        lines!(ax1, x, y, color = t, colormap = :viridis, linewidth = 0.5)
        ax2 = Axis(fig2[1,1], ylabel = L"S_N", title = "Wood's 5-box model'")
        lines!(ax2, t, x, color = :red, linewidth = 1.0)
        ax3 = Axis(fig2[2,1], xlabel = L"t", ylabel = L"S_T")
        lines!(ax3, t, y, color = :blue, linewidth = 1.0)

        # Solve the stationary SDE
        ensemble = evolve(f, η, Λ_null, z0, stepsize=δt, steps=Nt)

        # Extract the solutions
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]
        H = ensemble.parameter
        t = ensemble.time

        # Plot the timeseries and phase space trajectory
        ax7 = Axis(fig3[1,1], ylabel = L"S_N", title = "Wood's 5-box model")
        lines!(ax7, t, x, color = :red, linewidth = 1.0)
        ax8 = Axis(fig3[2,1], xlabel = L"t", ylabel = L"S_T")
        lines!(ax8, t, y, color = :blue, linewidth = 1.0)

        # Compute and plot the histogram of the solution
        bins, pdf = fit_distribution(x, n_bins=n_bins)
        ax9 = Axis(fig3[1,2], xlabel = L"S_N", ylabel = "density")
        barplot!(ax9, bins, pdf, color = :red, strokecolor = :black, strokewidth = 1.0)
        bins, pdf = fit_distribution(y, n_bins=n_bins)
        ax10 = Axis(fig3[2,2], xlabel = L"S_T", ylabel = "density")
        barplot!(ax10, bins, pdf, color = :blue, strokecolor = :black, strokewidth = 1.0)

        #=
        V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)
        D = (std(x)^2)/2
        ρ(x, c) = exp(-V(x, c)/D)

        xs = +(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) - C[2])
        xu = -(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) + C[2])
        I = (-Inf,Inf)
        if xs > xu
                I = (xu, +Inf) 
        else
                I = (-Inf, xu) 
        end
        integral = IntegralProblem(ρ, (minimum(bins), maximum(bins)), C)
        quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
        p(x, c) = ρ(x, c)/(quadrature.u)

        domain = LinRange(minimum(bins), maximum(bins), 1000)
        lines!(ax, domain, [p(s, C) for s in domain], color = :goldenrod2, linewidth = 3.0)
        =#

        #-------------------#
        #    Saddle-node    #
        #-------------------#
        
        # Solve the ensemble slow-fast SDE
        equilibria = get_equilibria(F1, F2, μ0)
        z0 = [equilibria.stable[1]; μ0]
        ensemble = evolve(F, ζ, Λ, z0, stepsize=δt, steps=Nt)

        # Extract the solutions
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]
        M = ensemble.parameter
        t = ensemble.time

        # Plot the timeseries and phase space trajectory
        ax4 = Axis(fig1[1,2], xlabel = L"x", ylabel = L"y", title = "Saddle-node normal form")
        lines!(ax4, x, y, color = t, colormap = :viridis, linewidth = 0.5)
        ax5 = Axis(fig2[1,2], ylabel = L"x", title = "Saddle-node normal form")
        lines!(ax5, t, x, color = :red, linewidth = 1.0)
        ax6 = Axis(fig2[2,2], xlabel = L"t", ylabel = L"y")
        lines!(ax6, t, y, color = :blue, linewidth = 1.0)

        # Solve the stationary SDE
        ensemble = evolve(F, ζ, Λ_null, z0, stepsize=δt, steps=Nt)

        # Extract the solutions
        x = (ensemble.state[1])[:, 1]
        y = (ensemble.state[1])[:, 2]
        M = ensemble.parameter
        t = ensemble.time

        # Plot the timeseries and phase space trajectory
        ax11 = Axis(fig4[1,1], ylabel = L"x", title = "Saddle-node normal form")
        lines!(ax11, t, x, color = :red, linewidth = 1.0)
        ax12 = Axis(fig4[2,1], xlabel = L"t", ylabel = L"y")
        lines!(ax12, t, y, color = :blue, linewidth = 1.0)

        # Compute and plot the histogram of the solution
        bins, pdf = fit_distribution(x, n_bins=n_bins)
        ax13 = Axis(fig4[1,2], xlabel = L"S_N", ylabel = "density")
        barplot!(ax13, bins, pdf, color = :red, strokecolor = :black, strokewidth = 1.0)
        bins, pdf = fit_distribution(y, n_bins=n_bins)
        ax14 = Axis(fig4[2,2], xlabel = L"S_T", ylabel = "density")
        barplot!(ax14, bins, pdf, color = :blue, strokecolor = :black, strokewidth = 1.0)

        # Export figures
        savefig("comparison/phase_space.png", fig1)
        savefig("comparison/timeseries.png", fig2)
        savefig("comparison/OUP_Wood.png", fig3)
        savefig("comparison/OUP_SN.png", fig4)
end

# Execute the main
main()
