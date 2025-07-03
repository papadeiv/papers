include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timeseries data 
solution = readin("../data/slide_05/solution.csv")
t = solution[:,1]
u = solution[:,3]

# Import the ews data
data = readin("../data/slide_05/ews_variance.csv")
t_var = data[:,1]
u_var = data[:,2] 

# Loop over the realisations of the variance timeseries
printstyled("Generating figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:length(t_var)
        ##############
        # Timeseries #
        ##############

        # Create and customise the timeseries axis 
        fig, ax1 = mkfig(size = [1200,800],
                         bg_out = colorant"rgb(147,237,213)",
                         border_color = RGBAf(0,0.236,0.236),
                         limits = ((t[1], t[end]), (-1, 3)),
                         box_position = [1,1],
                         lab = [L"\textbf{time}", L"\textbf{state}"],
                         lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                         toggle_lab = [false, true],
                         toggle_ticks = [false,false],
                         toggle_ticks_lab = [false,false],
                        )
        # Plot the state timeseries 
        lines!(ax1, LinRange(0.0, t[end], length(t)), u, linewidth = 2, color = (:teal, 1.0))
        # Plot the graybox reppresenting the sliding window
        poly!(ax1, Point2f[(t[n], -1), (t_var[n], -1), (t_var[n], 3), (t[n], 3)], color = (colorant"rgb(147,237,213)", 0.35), strokecolor = :grey, strokewidth = 0.05)
 
        ############
        # Variance #
        ############

        # Create and customise the variance axis 
        Fig, ax2 = mkfig(fig = fig,
                         box_position = [2,1],
                         bg_out = colorant"rgb(147,237,213)",
                         border_color = RGBAf(0,0.236,0.236),
                         limits = ((t[1], t[end]), (-0.0005, 0.005)),
                         lab = [L"\textbf{time}", L"\textbf{variance}"],
                         lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                         toggle_ticks = [false,false],
                         toggle_ticks_lab = [false,false],
                        )
        # Plot the variance timeseries up until the current timestep 
        lines!(ax2, t_var[1:n], u_var[1:n], linewidth = 2, color = (:teal, 1.0))
        # Plot the graybox reppresenting the sliding window
        poly!(ax2, Point2f[(t[n], -0.0005), (t_var[n], -0.0005), (t_var[n], 0.02), (t[n], 0.02)], color = (colorant"rgb(147,237,213)", 0.35), strokecolor = :grey, strokewidth = 0.05)
        # Plot the value of the variance at the current timestep
        scatter!(ax2, t_var[n], u_var[n], color = :teal, strokecolor = :black, strokewidth = 1.5, markersize = 30)

        # Export the figure 
        save("../fig/slide_05/$n.png", fig)
end
