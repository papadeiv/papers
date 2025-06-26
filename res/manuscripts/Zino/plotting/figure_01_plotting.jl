include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

 # Import the initial conditions 
x0 = readin("../data/figure_01/x0.csv")
Nx = length(x0)

# Number of parameter realisations
Nμ = 3 

# Create and customise the figure 
fig, ax = mkfig(size = [2400,900],
                pad = (30,60,10,5), # Order is: left, right, bottom, top 
                bg_out = :white,
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false]
               )

# Loop over the parameter values
printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
         # Import the solution data at current timestep from csv 
        solution = readin("../data/figure_01/$n.csv")
        t = solution[:,1]
        X = solution[:,2:end]

        # Create and customise the figure 
        if n==1
                Fig, ax = mkfig(fig=fig,
                                title = L"\textbf{Coordination game}",
                                toggle_title = true,
                                title_size = 50,
                                box_position = [1,n],
                                limits = ((t[1], t[end]), (-0.01, 1.01)),
                                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                                lab_pad = [-60.0,-60.0],
                                x_ticks = [t[1],t[end]],
                                y_ticks = [0,1],
                                ticks_lab_xpos = [:right,:top],
                                ticks_lab_trunc = [0,1],
                               )
                # Define vector of indices for plotting purposes
                idx = collect(UnitRange(1::Int64, Nx))
                randomised = collect(UnitRange(1::Int64, Nx))
                StatsBase.direct_sample!(idx, randomised)
                # Loop over the different ICs
                for m in 1:Nx
                        # Extract the solution
                        x = X[:,m]
                        # Plot the solution 
                        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [:teal, :brown2], colorrange = (1,Nx))
                end
        elseif n==2
                Fig, ax = mkfig(fig=fig,
                                title = L"\textbf{Dominant strategy}",
                                toggle_title = true,
                                title_size = 50,
                                box_position = [1,n],
                                limits = ((t[1], t[end]), (-0.01, 1.01)),
                                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                                toggle_lab = [true,false],
                                lab_pad = [-60.0,-60.0],
                                x_ticks = [t[1],t[end]],
                                y_ticks = [0,1],
                                ticks_lab_xpos = [:center,:top],
                                toggle_ticks_lab = [true,false],
                                ticks_lab_trunc = [0,1],
                               )
                # Define vector of indices for plotting purposes
                idx = collect(UnitRange(1::Int64, Nx))
                randomised = collect(UnitRange(1::Int64, Nx))
                StatsBase.direct_sample!(idx, randomised)
                # Loop over the different ICs
                for m in 1:Nx
                        # Extract the solution
                        x = X[:,m]
                        # Plot the solution 
                        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [:teal, :brown2], colorrange = (1,Nx))
                end
        else
                Fig, ax = mkfig(fig=fig,
                                title = L"\textbf{Anti-coordination}",
                                toggle_title = true,
                                title_size = 50,
                                box_position = [1,n],
                                limits = ((t[1], t[end]), (-0.01, 1.01)),
                                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                                toggle_lab = [true,false],
                                lab_pad = [-60.0,-60.0],
                                x_ticks = [t[1],t[end]],
                                y_ticks = [0,1],
                                ticks_lab_xpos = [:left,:top],
                                toggle_ticks_lab = [true,false],
                                ticks_lab_trunc = [0,1],
                               )
                # Define vector of indices for plotting purposes
                idx = collect(UnitRange(1::Int64, Nx))
                randomised = collect(UnitRange(1::Int64, Nx))
                StatsBase.direct_sample!(idx, randomised)
                # Loop over the different ICs
                for m in 1:Nx
                        # Extract the solution
                        x = X[:,m]
                        # Plot the solution 
                        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [:teal, :brown2], colorrange = (1,Nx))
                end
        end
end

# Export the figure 
save("../fig/figure_01.png", fig)
