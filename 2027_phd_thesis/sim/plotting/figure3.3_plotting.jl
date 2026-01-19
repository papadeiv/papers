include("../../../inc/PlottingTools.jl")
include("../../../inc/IO.jl")

# Import the timseries data 
data = readin("./data/timeseries.csv")
time = data[:,1]
state = data[:,2]

# Define quantities for the plots
trueVar = 0.005
∂N = 0.35

# Create the figure and the axis for the plot 
fig, ax1 = mkfig(size = [2000,1000], lab = [L"\textbf{time}",L"\mathbf{x}"], box_position = [1,1:3])
nullfig, ax2 = mkfig(fig = fig, 
                     box_position = [1,4], 
                     lab = [L"\textbf{density}", L"\mathbf{x}"], 
                     toggle_lab = [true, false],
                     x_ticks = [0.00, 0.06],
                     y_ticks = [-∂N, ∂N], 
                     toggle_ticks_lab = [true, false], 
                     ticks_lab_trunc = [2,0]
                    )

# Setup the ticks for the axis
set_ticks(ax1, time, state)

# Plot the timeseries
lines!(ax1, time, state, color = :black, linewidth = 3)
hist!(ax2, state, bins = 50, normalization = :probability, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)

# Plot the critical manifold 
lines!(ax1, [time[1], time[end]], [0.0, 0.0], color = (:grey,0.5), linewidth = 50)
# Plot the boundaries of the neighbour of the critical manifold
lines!(ax1, [time[1], time[end]], [∂N, ∂N], color = :blue, linewidth = 3, linestyle = :dash)
text!(ax1, 0.0, +0.4, text=L"\mathbf{\partial N^{-}(\rho;C_{\epsilon})}", color = :blue, align = (:left, :center), fontsize = 30)
lines!(ax1, [time[1], time[end]], [-∂N, -∂N], color = :blue, linewidth = 3, linestyle = :dash)
text!(ax1, 0.0, -0.4, text=L"\mathbf{\partial N^{+}(\rho;C_{\epsilon})}", color = :blue, align = (:left, :center), fontsize = 30)

# Export the figure
save("./fig/timeseries.png", fig)
