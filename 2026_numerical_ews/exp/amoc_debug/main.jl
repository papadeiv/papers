"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")

# Define the main algorithm
function main()
        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/Atlantic_250-500m_year_0001-2200.nc")
        time = data["time"][:]
        latitude = data["lat"][:]
        temperature = data["TEMP"]
        salinity = data["SALT"]

        N = 12::Integer
        t = copy(time[1:idx])
        T = copy(temperature[N,1:idx]) 
        S = copy(salinity[N,1:idx]) 

        # Define the observable (flow)
        ϕ = T - S
        detrended_solution = detrend(ϕ, alg = "emd", n_modes=1)
        trend = detrended_solution.trend
        residual = detrended_solution.residuals

        ensemble = preprocess_solution(t, residual, window_size)
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Loop over the ensemble
        ews = Vector{Float64}(undef, Ne)
        for n in 1:Ne
                # Generate figure and axis
                fig = Figure(size = (1200,800))
                ax1 = Axis(fig[1,1:2], xlabel = L"\textbf{years}", ylabel = L"\textbf{flow (250-500m)}", limits = (time[1], time[end], -32,-20))
                ax2 = Axis(fig[2,1:2], xlabel = L"\textbf{years}", ylabel = L"\textbf{early-warning signal}", limits = (time[1], time[end], -0.1, 1.1))
                ax3 = Axis(fig[1,3], xlabel = L"\textbf{residuals}", ylabel = L"\textbf{density}", limits = (-0.2, 0.2, nothing, 15))
                ax4 = Axis(fig[2,3], xlabel = L"\mathbf{V(x)}", ylabel = L"\mathbf{x}", limits = (-5,5,-5,5))

                # Extract windowed subseries
                tw = copy(ensemble.timesteps[n])
                ϕw = copy(ensemble.trajectories[n]) 

                # Assemble histogram of the windowed subseries
                σ = std(ϕw)
                n_bins = convert(Int64, ceil(abs(maximum(ϕw)-minimum(ϕw))/(3.49*σ*(length(ϕw))^(-1.0/3.0))))
                bins, pdf = fit_distribution(ϕw, n_bins=n_bins)

                println("$n) std = $σ, n_bins = $(n_bins)")

                # Solve the nonlinear least-squares problem and compute the early-warning signal
                solution = fit_potential(ϕw, transformation=[0.0,0.0,16.0], n_bins=n_bins).fit
                ews[n] = analyse(solution)

                # Plot and export the results
                interval = LinRange(minimum(bins), maximum(bins), 1000)
                domain = LinRange(-5,5,1000)
                lines!(ax1, t, ϕ, color = (:red,0.35), linewidth = 1.0)
                lines!(ax1, time[(idx+1):end], temperature[N,(idx+1):end] - salinity[N,(idx+1):end], color = (:black,0.35), linewidth = 1.0)
                lines!(ax1, t[n:(Nt+n-1)], ϕ[n:(Nt+n-1)], color = :red, linewidth = 1.0)
                poly!(ax1, Point2f[(tw[1], -32), 
                                   (tw[end], -32),
                                   (tw[end], -20),
                                   (tw[1], -20),
                                  ], color = (CtpGray, 0.1), strokecolor = :black, strokewidth = 1.0)
                lines!(ax1, [time[idx], time[idx]], [-32, -20], color = :red, linestyle = :dash, linewidth = 2.0)
                lines!(ax2, t[(Nt-1):(Nt+n-2)], ews[1:n], color = :red, linewidth = 4.0)
                lines!(ax2, [time[idx], time[idx]], [-32, -20], color = :red, linestyle = :dash, linewidth = 2.0)
                scatter!(ax2, t[(Nt+n-2)], ews[n], color = :red, markersize = 20.0, strokewidth = 1.0, strokecolor = :black)
                barplot!(ax3, bins, pdf, color = pdf, colormap = [:white,(:red,0.5)], strokecolor = :black, strokewidth = 1.0)
                I = (-Inf,Inf)
                if xs(solution) > xu(solution)
                        I = (xu(solution), +Inf) 
                else
                        I = (-Inf, xu(solution)) 
                end
                ρ(x, c) = exp(-2*V(x, c)/σ^2)
                integral = IntegralProblem(ρ, I, solution)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                Z = 1/(quadrature.u)
                p(x, c) = Z*ρ(x, c)
                lines!(ax3, interval, [p(x, solution) for x in interval], color = :red, linewidth = 4.0)
                lines!(ax4, domain, [V(x, solution) for x in domain], color = :red, linewidth = 2.0)
                savefig("amoc/$(latitude[N])/$n.png", fig)
        end
end

# Execute the main
main()
