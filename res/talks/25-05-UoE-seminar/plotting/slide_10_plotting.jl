include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Import the data from csv
parameters = readin("../data/slides_10_to_13/parameters.csv")
x_stb = parameters[:,3]
solution = readin("../data/slides_10_to_13/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(t)

# Define the width of the sliding window (as a relative length of the timeseries)
width = 0.150::Float64

# Assemble the sliding window
window = get_window_parameters(Nt, width)
Ns = window[2]

# Loop over the window's strides 
printstyled("Generating the figures for the histograms using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Ns
        ##############
        # Timeseries #
        ##############

        # Import the data
        u_det = readin("../data/slide_10/residuals/$n.csv")
        t_loc = u_det[:,1]
        μ_loc = u_det[:,2]
        u_loc = u_det[:,3]
        mean_trend = u_det[:,4]
        mean_residuals = u_det[:,5]
        lin_trend = u_det[:,6]
        lin_residuals = u_det[:,7]

        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64
        # Derive the drift of the QSE in the sliding window
        x_loc = x_stb[n:n_end]

        # Get range of the timeseries for plotting
        x_inf = minimum(u_loc)
        x_sup = maximum(u_loc)
        mean_res_inf = minimum(mean_residuals)
        mean_res_sup = maximum(mean_residuals)
        lin_res_inf = minimum(lin_residuals)
        lin_res_sup = maximum(lin_residuals)

        # Create and customise the non-stationary timeseries figure 
        fig, ax = mkfig(size = [1500,1200],
                        box_position = [1,1:2],
                        limits = ((t_loc[1],t_loc[end]), (x_inf,x_sup)),
                        toggle_lab = [false,true],
                        lab = [L"\textbf{time}", L"\mathbf{x(t)}"],
                        lab_pad = [-50.0,-50.0],
                        x_ticks = [t_loc[1],t_loc[end]],
                        y_ticks = [x_inf,x_sup],
                        toggle_ticks_lab = [false,true],
                        ticks_lab_trunc = [0,1]
                       )
        # Plot the non-stationary timeseries
        lines!(ax, t_loc, u_loc, linewidth = 2.5, color = :teal)
        # Plot the drift of the QSE
        lines!(ax, t_loc, x_loc, linewidth = 5, color = :black)
        # Plot the mean trend 
        lines!(ax, t_loc, mean_trend, linewidth = 4.5, color = colorant"rgb(207,19,92)")
        # Plot the linear trend 
        lines!(ax, t_loc, lin_trend, linewidth = 4.5, color = colorant"rgb(38,106,176)")

        # Create and customise the mean detrended timeseries figure 
        fig, ax = mkfig(fig = fig,
                        box_position = [2,1:2],
                        limits = ((t_loc[1],t_loc[end]), (mean_res_inf,mean_res_sup)),
                        toggle_lab = [false,true],
                        lab = [L"\textbf{time}", L"\textbf{residuals}"],
                        lab_pad = [-50.0,-65.0],
                        x_ticks = [t_loc[1],t_loc[end]],
                        y_ticks = [mean_res_inf,mean_res_sup],
                        toggle_ticks_lab = [false,true],
                        ticks_lab_trunc = [0,1]
                       )
        # Plot the mean detrended timeseries
        lines!(ax, t_loc, mean_residuals, linewidth = 2.5, color = colorant"rgb(207,19,92)")
 
        # Create and customise the linear detrended timeseries figure 
        fig, ax = mkfig(fig = fig,
                        box_position = [3,1:2],
                        limits = ((t_loc[1],t_loc[end]), (lin_res_inf,lin_res_sup)),
                        lab = [L"\textbf{time}", L"\textbf{residuals}"],
                        lab_pad = [-50.0,-65.0],
                        x_ticks = [t_loc[1],t_loc[end]],
                        y_ticks = [lin_res_inf,lin_res_sup],

                        ticks_lab_xpos = [:right,:top],
                        ticks_lab_trunc = [0,1]
                       )
        # Plot the linear detrended timeseries
        lines!(ax, t_loc, lin_residuals, linewidth = 2.5, color = colorant"rgb(38,106,176)")
 
        #########################
        # Approximate histogram #
        #########################

        # Import the data
        mean_hist = readin("../data/slide_10/histograms/SN_mean/$n.csv")
        mean_bins = mean_hist[:,1]
        mean_pdf = mean_hist[:,2]
        lin_hist = readin("../data/slide_10/histograms/SN_lin/$n.csv")
        lin_bins = lin_hist[:,1]
        lin_pdf = lin_hist[:,2]

        # Create and customise the mean detrended histogram figure 
        fig, ax = mkfig(fig = fig,
                        box_position = [2,3],
                        limits = ((-0.05,6), (mean_res_inf,mean_res_sup)),
                        toggle_lab = [false,false],
                        lab = [L"\mathbf{\rho}",L"\mathbf{x}"],
                        lab_pad = [-50.0,-50.0],
                        x_ticks = [0,6],
                        y_ticks = [mean_res_inf,mean_res_sup],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the mean detrended histogram 
        lines!(ax, mean_pdf, mean_bins, linewidth = 4.5, color = colorant"rgb(207,19,92)")

        # Create and customise the linear detrended histogram figure 
        fig, ax = mkfig(fig = fig,
                        box_position = [3,3],
                        limits = ((-0.05,6), (lin_res_inf,lin_res_sup)),
                        toggle_lab = [true,false],
                        lab = [L"\mathbf{\rho}",L"\mathbf{x}"],
                        lab_pad = [-60.0,-50.0],
                        x_ticks = [0,6],
                        y_ticks = [lin_res_inf,lin_res_sup],
                        ticks_lab_xpos = [:left,:top],
                        toggle_ticks_lab = [true,false],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the linear detrended histogram 
        lines!(ax, lin_pdf, lin_bins, linewidth = 4.5, color = colorant"rgb(38,106,176)")

        # Export the figure
        save("../fig/slide_10.$n.png", fig)
end
