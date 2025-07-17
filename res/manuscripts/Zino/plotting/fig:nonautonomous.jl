include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

#=
# Import the initial conditions 
x0 = readin("../data/fig:nonautonomous/x0.csv")
Nx = length(x0)

# Import the parameter shift rates 
rates = readin("../data/fig:nonautonomous/rates.csv")
Nε = length(rates)

# Import the shifted unstable equilibria 
x3 = readin("../data/fig:nonautonomous/x3.csv")
Nt = length(x3[:,1])

# Import the R-tipping statistics 
tip_cnt = convert(AbstractArray{Int64}, readin("../data/fig:nonautonomous/tip_cnt.csv"))
tip_idx = convert(AbstractArray{Int64}, readin("../data/fig:nonautonomous/tip_idx.csv"))

# Loop over the rates
printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε
        # Get the current rate
        ε = rates[n]

        # Import the solutions at the current parameter shift rate 
        solutions = readin("../data/fig:nonautonomous/solutions/$n.csv")
        local t = solutions[:,1]
        local λ = solutions[:,2]
        local X = solutions[:,3:end]

        #######################
        #   Bifurcation set   #
        #######################

        # Define the speed of the parameter shift
        k = 0.5
        δ = 0.0
        Λ(t) = k*ε*(sech(ε*t + δ))^2

        # Compute the time realisations of ɑ and β
        β = λ
        ɑ = 1.0 .- β
        
        # Define limits of the set
        ɑ_min = -(0.1::Float64)
        ɑ_max = 2.0::Float64
        β_min = -(0.1::Float64)
        β_max = 2.0::Float64

        # Create and customise the figure 
        fig, ax = mkfig(size = [2100,700],
                        bg_out = :white,
                        pad = (20,25,10,20), # Order is: left, right, bottom, top 
                        box_position = [1,3],
                        limits = ((ɑ_min,ɑ_max), (β_min,β_max)),
                        lab = [L"\mathbf{\alpha:=a-c}",L"\mathbf{\beta:=d-b}"],
                        lab_pad = [-60.0,-30.0],
                        x_ticks = [0.0,ɑ_max],
                        y_ticks = [0.0,β_max],
                        ticks_lab_trunc = [0,0],
                       )
        # Plot the separatrices of the bifurcation set 
        lines!(ax, [0,0], [β_min,β_max], linewidth = 6, color = :black)
        lines!(ax, [ɑ_min,ɑ_max], [0,0], linewidth = 6, color = :black)
        # Plot the realisations in the bifurcation set
        lines!(ax, ɑ, β, linewidth = 20, color = [Λ(T) for T in t], colormap = :viridis)
        # Plot the start and endpoint of the trajectory in the bifurcation set
        scatter!(ax, ɑ[1], β[1], color = :brown2, markersize = 40, strokecolor = :black, strokewidth = 2)
        scatter!(ax, ɑ[end], β[end], color = :darkgreen, markersize = 40, strokecolor = :black, strokewidth = 2)
        # Plot the colorbar for all the rates
        Colorbar(fig[1,4],
                 size = 30, 
                 spinewidth = 5.0,
                 limits = (Λ(t[1]), Λ(0.0)), 
                 colormap = (:viridis, 1.0),
                 ticks = [Λ(t[1]), Λ(0.0)],
                 ticksize = 22,
                 tickwidth = 5.0,
                 tickformat = "{:.$(2)f}",
                 ticklabelsize = 50,
                 label = L"\mathbf{\dot{\lambda}(t)}",
                 labelpadding = -80.0,
                 labelsize = 50
                )

        #########################
        #   Coordination game   #
        #########################

        # Extract tipped indices of the set of solutions at the current shift rate
        indices = tip_idx[1:tip_cnt[n],n]
        
        # Create and customise the figure 
        fig, ax = mkfig(fig = fig,
                        box_position = [1,1:2],
                        limits = ((t[1], t[end]), (-0.05, 1.05)),
                        lab = [L"\mathbf{t}",L"\mathbf{x}"],
                        lab_pad = [-60.0,-30.0],
                        x_ticks = [t[1],t[end]],
                        y_ticks = [0,1],
                        ticks_lab_trunc = [0,0],
                       )
        # Loop over the different ICs
        for m in 1:Nx
                # Extract the solution
                x = X[:,m]
               # Check if solution at the current index has tipped
                if any(s -> s == m, indices) 
                        # Plot the solution in red 
                        lines!(ax, t, x, linewidth = 4, color = (:brown2,0.35))
                else
                        # Plot the solution in green 
                        lines!(ax, t, x, linewidth = 4, color = (:darkgreen,0.35))
                end
        end
        # Plot the equilibria
        lines!(ax, [t[1], t[end]], [0.0, 0.0], linewidth = 6, color = :black)
        lines!(ax, [t[1], t[end]], [1.0, 1.0], linewidth = 6, color = :black)
        lines!(ax, t, x3[:,n], linewidth = 6, color = :black, linestyle = :dash)

        # Export the figure 
        save("../fig/fig:nonautonomous/solutions/$n.png", fig)
end
=#

# Import an example solution
solutions = readin("../data/fig:nonautonomous/solutions/1.csv")
t = solutions[:,1]

# Create and customise the figure for the parameter shift
fig, ax1 = mkfig(size = [1400,700],
                 bg_out = :white,
                 pad = (15,35,10,20), # Order is: left, right, bottom, top 
                 box_position = [1,1],
                 limits = ((t[1],t[end]), (-0.05,1.05)),
                 lab = [L"\mathbf{t}",L"\mathbf{\Lambda(t)}"],
                 toggle_lab = [false,true],
                 lab_pad = [-60.0,-30.0],
                 x_ticks = [t[1],t[end]],
                 y_ticks = [0.0,1.0],
                 toggle_ticks_lab = [false,true],
                 ticks_lab_trunc = [0,0],
                )
 
# Create and customise the figure for the time derivative of the shift 
fig, ax2 = mkfig(fig = fig,
                 box_position = [2,1],
                 limits = ((t[1],t[end]), (-0.01,0.2)),
                 lab = [L"\mathbf{t}",L"\mathbf{\dot{\Lambda}(t)}"],
                 lab_pad = [-60.0,-30.0],
                 x_ticks = [t[1],t[end]],
                 y_ticks = [0.0,0.2],
                 ticks_lab_trunc = [0,0],
                )

# Loop over the rates
for n in 1:Nε
        # Get the current rate
        ε = rates[n]

        # Define the parameter shift and its time derivative
        Λ(t) = 0.5*(tanh(ε*t) + 1.0)
        Λt(t) = 0.5*ε*(sech(ε*t))^2

        # Plot both timeseries
        lines!(ax1, t, [Λ(T) for T in t], linewidth = 5, color = ε, colormap = (:viridis, 0.5), colorrange = (rates[1],rates[end]))
        lines!(ax2, t, [Λt(T) for T in t], linewidth = 5, color = ε, colormap = (:viridis, 0.5), colorrange = (rates[1],rates[end]))
end

# Plot the colorbar for all the rates
Colorbar(fig[1:2,2],
         size = 30, 
         spinewidth = 5.0,
         limits = (rates[1], rates[end]), 
         colormap = (:viridis, 1.0),
         ticks = [rates[1], rates[end]],
         ticksize = 22,
         tickwidth = 5.0,
         tickformat = "{:.$(2)f}",
         ticklabelsize = 50,
         label = L"\mathbf{\varepsilon}",
         labelpadding = -80.0,
         labelsize = 50
)

# Export the figure 
save("../fig/fig:nonautonomous/fig:sigmoid_example.png", fig)
