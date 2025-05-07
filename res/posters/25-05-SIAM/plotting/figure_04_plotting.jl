include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define window values to plot
init = 301 
fin = 2419 

# Import the timeseries data 
Time = readin("../data/figure_04/time.csv")
time = Time[init:fin]
U = readin("../data/figure_04/u.csv")
u = U[init:fin]

# Import the ews data
data = readin("../data/figure_04/ews_variance.csv")
T_var = data[:,1]
t_var = T_var[1:(fin-init)]
U_var = data[:,2] 
u_var = U_var[1:(fin-init)]

# Create and customise the timeseries axis 
fig, ax = mkfig(size = [1600,800],
                border = 10.0,
                limits = ((time[1], time[end]), (1.75, 2.5)),
                box_position = [1,1],
                lab = [L"\textbf{time}", L"\mathbf{x_t}"],
                toggle_lab = [true, true],
                lab_color = [:black,:darkgoldenrod2],
                lab_size = [70,70],
                lab_pad = [10.0,-60.0],
                x_ticks = [time[1],time[end]],
                y_ticks = [1.75, 2.5],
                toggle_ticks = [false,true],
                toggle_ticks_lab = [false,true],
                ticks_lab_size = [70,70],
                ticks_lab_trunc = [0,2]
)
# Plot the state timeseries 
lines!(ax, LinRange(0.0, time[end], length(time)), u, linewidth = 6, color = (:darkgoldenrod2, 1.0))
 
# Mirror the axis for the variance
ax = mirror_axis(fig, [time[1],time[end]];
                 color = :teal,
                 y_lab = L"\textbf{var(}\mathbf{x_t}\textbf{)}",
                )
ax.limits = ((time[1], time[end]), (0, 0.02))
ax.yticks = [0, 0.02] 
ax.ylabelsize = 70.0
ax.yticklabelsize = 70
ax.ytickwidth = 10.0

# Plot the variance timeseries up until the current timestep 
lines!(ax, t_var, u_var, linewidth = 8, color = (:teal, 1.0))

# Export the figure 
save("../fig/fig4.png", fig)
