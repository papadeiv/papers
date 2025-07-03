include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timeseries data 
solution = readin("../data/slide_03/solution.csv")
t = solution[:,1]
u = solution[:,3]

# Import the ews data
data = readin("../data/slide_03/ews.csv")
t_var = data[:,1]
u_var = data[:,2] 
u_skw = data[:,3] 
Ns = length(t_var)

# Define a list of plotting indices
idx = collect(StepRange(1, Int64(20), Ns)) 

# Loop over the realisations of the variance timeseries
printstyled("Generating figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for m in 1:length(idx)
        # Get the current plotting index
        n = idx[m] 
 
        ##############
        # Timeseries #
        ##############

        # Create and customise the timeseries axis 
        fig, ax1 = mkfig(size = [900,600],
                         bg_out = colorant"rgb(239,241,245)",
                         border_color = colorant"rgb(76,79,105)",
                         limits = ((t[1], t[end]), (-1, 3)),
                         box_position = [1,1],
                         lab = [L"\textbf{time}", L"\textbf{state}"],
                         lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                         toggle_lab = [false, true],
                         toggle_ticks = [false,false],
                         toggle_ticks_lab = [false,false],
                        )
        # Plot the state timeseries 
        lines!(ax1, LinRange(0.0, t[end], length(t)), u, linewidth = 2, color = (colorant"rgb(23,146,153)", 1.0))
        # Plot the graybox reppresenting the sliding window
        poly!(ax1, Point2f[(t[n], -1), (t_var[n], -1), (t_var[n], 3), (t[n], 3)], color = (colorant"rgb(156,160,176)", 0.35), strokecolor = :grey, strokewidth = 0.05)
 
        ###########
        #   EWS   #
        ###########

        # Create and customise the variance axis 
        Fig, ax2 = mkfig(fig = fig,
                         box_position = [2,1],
                         bg_out = colorant"rgb(239,241,245)",
                         border_color = colorant"rgb(76,79,105)",
                         limits = ((t[1], t[end]), (0.0005, 0.004)),
                         lab = [L"\textbf{time}", L"\textbf{variance}"],
                         lab_color = [colorant"rgb(76,79,105)", colorant"rgb(136,57,239)"],
                         toggle_ticks = [false,false],
                         toggle_ticks_lab = [false,false],
                        )
        # Plot the graybox reppresenting the sliding window
        poly!(ax2, Point2f[(t[n], -0.0005), (t_var[n], -0.0005), (t_var[n], 0.02), (t[n], 0.02)], color = (colorant"rgb(156,160,176)", 0.35), strokecolor = :grey, strokewidth = 0.05)
        # Plot the variance timeseries up until the current timestep 
        lines!(ax2, t_var[1:n], u_var[1:n], linewidth = 3, color = (colorant"rgb(136,57,239)", 1.0))
        # Plot the value of the variance at the current timestep
        scatter!(ax2, t_var[n], u_var[n], color = colorant"rgb(136,57,239)", strokecolor = :black, strokewidth = 1.5, markersize = 30)

        # Create and customise the skewness mirrored axis
        ax3 = mirror_axis(fig, [t[1], t[end]], [-0.75, 0.5],
                          box_position = [2,1],
                          color = colorant"rgb(230,69,83)",
                          y_lab = L"\textbf{skewness}",
                          toggle_ticks = false 
                         )
        # Plot the skewness timeseries up until the current timestep 
        lines!(ax3, t_var[1:n], u_skw[1:n], linewidth = 3, color = (colorant"rgb(230,69,83)", 1.0))
        # Plot the value of the skewness at the current timestep
        scatter!(ax3, t_var[n], u_skw[n], color = colorant"rgb(230,69,83)", strokecolor = :black, strokewidth = 1.5, markersize = 30)

        # Export the figure 
        save("../fig/slide_03/$m.png", fig)
end
