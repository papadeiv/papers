include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timeseries data 
time = readin("../data/slide_07/time.csv")
u = readin("../data/slide_07/u.csv")

# Import the ews data
data = readin("../data/slide_07/ews_variance.csv")
t_var = data[:,1]
u_var = data[:,2] 

# Loop over the realisations of the variance timeseries
printstyled("Generating figures\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(t_var)
        ##############
        # Timeseries #
        ##############

        # Create and customise the timeseries axis 
        fig, ax1 = mkfig(size = [1200,800],
                bg_out = "#eeeeeeff",
                limits = ((time[1], time[end]), (-1, 3)),
                box_position = [1,1],
                lab = [L"\textbf{time}", L"\textbf{state}"],
                toggle_lab = [false, true],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
        # Plot the state timeseries 
        lines!(ax1, LinRange(0.0, time[end], length(time)), u, linewidth = 3.5, color = (:black, 1.0))
        # Plot the graybox reppresenting the sliding window
        poly!(ax1, Point2f[(time[n], -1), (t_var[n], -1), (t_var[n], 3), (time[n], 3)], color = (:grey, 0.35), strokecolor = :grey, strokewidth = 0.05)
 
        ############
        # Variance #
        ############

        # Create and customise the variance axis 
        Fig, ax2 = mkfig(fig = fig,
                             box_position = [2,1],
                             limits = ((time[1], time[end]), (-0.0005, 0.02)),
                             lab = [L"\textbf{time}", L"\textbf{variance}"],
                             toggle_ticks = [false,false],
                             toggle_ticks_lab = [false,false],
                    )
        # Plot the variance timeseries up until the current timestep 
        lines!(ax2, t_var[1:n], u_var[1:n], linewidth = 3.5, color = (:brown2, 1.0))
        # Plot the graybox reppresenting the sliding window
        poly!(ax2, Point2f[(time[n], -0.0005), (t_var[n], -0.0005), (t_var[n], 0.02), (time[n], 0.02)], color = (:grey, 0.35), strokecolor = :grey, strokewidth = 0.05)
        # Plot the value of the variance at the current timestep
        scatter!(ax2, t_var[n], u_var[n], color = :brown2, strokecolor = :black, strokewidth = 1.5, markersize = 30)

        # Export the figure 
        save("../fig/slide_07/$n.png", fig)
end
