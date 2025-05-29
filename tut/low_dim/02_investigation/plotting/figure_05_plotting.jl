include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timestamps and parameter values
time = readin("../data/figure_05/time.csv")
Nt = length(time)
parameter = readin("../data/figure_05/parameter.csv")
Nμ = length(parameter)
metrics = readin("../data/figure_05/metrics.csv")
ews = readin("../data/figure_05/variance.csv")

# Define the spatial domain
L = 2*π 
Nx = convert(Int64, 1e3)
dx = L/Nx
I = LinRange(-L, L, Nx)
r = 0.2

# Index of images you want to export
idx = [1, 459, 1000]

# Loop over the parameter sweep
printstyled("Generating the figures\n"; bold=true, underline=true, color=:green)
@showprogress for n in 1:length(idx)
        # Import the solution
        solution = readin("../data/figure_05/solutions/$(idx[n]).csv") 

        # Format the parameter value
        μ = trunc(parameter[idx[n]], digits=3)

        # Create and customise the timeseries axis 
        fig, ax = mkfig(size = [1000,1000],
                        border = 15.0,
                        title = L"\text{\mu=%$(μ)}",
                        toggle_title = true,
                        title_size = 120,
                        title_gap = 14.0,
                        limits = ((-r*L,r*L), (-1.5,1.5)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false]
                        )

        # Plot the steady-state
        lines!(ax, I, solution[:,end], linewidth = 8.0, color = (:teal, 1.0))

        # Export the figure 
        save("../fig/fig5.$n.png", fig)
end

# Create and customise the metrics axis 
fig, ax = mkfig(size = [1600,800],
                border = 10.0,
                limits = ((parameter[1], parameter[end]), (0,40)),
                lab = [L"\textbf{\mu}", L"\mathbf{\bar{u}}"],
                toggle_lab = [true, true],
                lab_size = [70,70],
                lab_color = [:black, :darkgoldenrod2],
                lab_pad = [-80.0, -80.0],
                x_ticks = [parameter[1],parameter[end]],
                y_ticks = [0, 40],
                toggle_ticks = [true,true],
                toggle_ticks_lab = [true,true],
                ticks_lab_size = [70,70],
                ticks_lab_trunc = [1,0]
)
# Plot the spatial mean 
lines!(ax, parameter, metrics[:,1], linewidth = 7.0, color = (:darkgoldenrod2, 1.0))

# Mirror the axis for the variance
ax = mirror_axis(fig, [time[1],time[end]];
                 color = :teal,
                 y_lab = L"\textbf{var(}\mathbf{\bar{u}}\textbf{)}",
                )
ax.limits = ((parameter[1], parameter[end]), (-0.0000025, 0.02))
ax.ylabelsize = 70
ax.ylabelpadding = -80.0
ax.yticks = [0, 0.02] 
ax.yticklabelsize = 70 
ax.ytickformat = "{:.$(2)f}"
ax.ytickwidth = 10.0

# Plot the variance timeseries up until the current timestep 
lines!(ax, ews[:,1], ews[:,2], linewidth = 7.0, color = (:teal, 1.0))

# Export the figure 
save("../fig/fig5.png", fig)
