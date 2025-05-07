include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timestamps and parameter values
time = readin("../data/figure_05/time.csv")
Nt = length(time)
solutions = readin("../data/figure_05/solutions.csv")
Ne = length(solutions[:,1])
escapes = readin("../data/figure_05/escapes.csv")
N_esc = length(escapes[:,1])
tipping = readin("../data/figure_05/tipping.csv")
N_tip = length(tipping[:,1])
T_dist = readin("../data/figure_05/escape_time_distribution.csv")
N_bins = length(T_dist[:,1])
mean_escaped = readin("../data/figure_05/mean_escaped.csv")
N_mean = length(mean_escaped)

# Define the plot limits for the equilibrium distribution histogram 
x_inf = -(1.125::Float64)
x_sup = 3.982::Float64
y_inf = -(2.00::Float64)
y_sup = 3.00::Float64 
domain = LinRange(x_inf, x_sup, 1000)

# Define the parameter value
μ = 0.000::Float64

# Define the potential function
V(x) = μ*x + x^2 - x^3 + (1/5)*x^4

# Create and customise the potential axis 
fig, ax = mkfig(size = [1600,1600],
                box_position = [1:2,1:3],
                limits = ((x_inf, x_sup), (y_inf, y_sup)),
                toggle_lab = [false,true],
                lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                lab_size = [40,40],
                lab_pad = [10.0, 10.0],
                x_ticks = [x_inf, x_sup],
                y_ticks = [y_inf, y_sup],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
                ticks_lab_trunc = [1,0]
)
# Plot the potential function 
lines!(ax, domain, [V(x) for x in domain], color = (:darkgoldenrod2,1.00), linewidth = 6)
# Plot the IC 
scatter!(ax, solutions[1,1], V(solutions[1,1]), color = :teal, markersize = 30, strokewidth = 2)
# Plot dashed line connecting the plot below
lines!(ax, [0.867,0.867], [y_inf, V(0.867)], color = (:brown2,1.00), linewidth = 5, linestyle = :dash)
# Plot the marker of the unstable equilibrium
scatter!(ax, 0.867, V(0.867), color = :brown2, marker = :star8, markersize = 35, strokewidth = 2)

# Create and customise the ensemble timeseries axis 
fig, ax = mkfig(fig = fig,
                box_position = [3:4,1:3],
                limits = ((x_inf,x_sup), (time[1], time[end])),
                toggle_lab = [true,true],
                lab = [L"\mathbf{x}", L"\textbf{time}"],
                lab_size = [40,40],
                lab_pad = [20.0, 10.0],
                ax_orientation = [false,true],
                x_ticks = [x_inf, x_sup],
                y_ticks = [time[1], time[250], time[500], time[end]],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
                ticks_lab_trunc = [0,2]
)
# Loop over the tipped ensemble
printstyled("Generating the ensemble timeseries subfigure\n"; bold=true, underline=true, color=:darkgreen)
@showprogress for e in 1:N_tip
        # Extract the trajectory index from the ensemble
        E_idx = convert(Int64, tipping[e,1])
        # Plot the timeseries 
        lines!(ax, solutions[E_idx,:], time, color = (:teal,0.025), linewidth = 2)
end
# Plot a mean escaped timeseries 
lines!(ax, solutions[mean_escaped[6],:], time, color = (:darkgoldenrod2,1.00), linewidth = 5)
# Plot dashed line connecting the plot above 
lines!(ax, [0.867,0.867], [0, escapes[1644,2]], color = (:brown2,1.00), linewidth = 5, linestyle = :dash)
# Plot dashed line connecting the plot on the right 
lines!(ax, [0.867,x_sup], [escapes[1644,2], escapes[1644,2]], color = (:brown2,1.00), linewidth = 5, linestyle = :dash)
# Mark its escape time
scatter!(ax, 0.867, escapes[1644,2], color = :brown2, marker = :star8, markersize = 35, strokewidth = 2)

# Create and customise the escape time axis 
fig, ax = mkfig(fig = fig,
                box_position = [3:4,4:5],
                limits = ((-0.01,0.7), (time[1], time[end])),
                toggle_lab = [true,false],
                lab = [L"\mathbf{p_{esc}(t)}", L"\textbf{time}"],
                lab_size = [40,40],
                lab_pad = [10.0, 10.0],
                ax_orientation = [false,true],
                flip_y = true,
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
)
# Plot the escape time histogram 
lines!(ax, T_dist[:,2], T_dist[:,1], color = (:teal,1.00), linewidth = 6)
lines!(ax, [0,0], [T_dist[end,1],time[end]], color = (:teal,1.00), linewidth = 6)
# Plot dashed line connecting the plot on the right 
lines!(ax, [0,0.305], [escapes[1644,2], escapes[1644,2]], color = (:brown2,1.00), linewidth = 5, linestyle = :dash)
# Mark the escape probability of the target timeseries
scatter!(ax,0.305, escapes[1644,2], color = :brown2, marker = :star8, markersize = 35, strokewidth = 2)

# Export the figure 
save("../fig/fig5.png", fig)
