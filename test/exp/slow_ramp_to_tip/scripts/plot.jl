"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"
CtpYellow = colorant"rgb(229,200,144)"
CtpWhite = colorant"rgb(198,208,245)"
CtpGray = colorant"rgb(98,104,128)"

#---------------------#
#    Tipping point    #
#---------------------#

# Figure for the timeseries 
fig1, ax1 = makefig(size = [1600,700],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    #title = L"\textbf{N_t = %$Nt,\;N_e = %$Ne}",
                    #toggle_title = true,
                    #limits = ((0, nothing), (-0.3, 0.3)),
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    #y_ticks = [-0.3,0.3],
                    ticks_lab_trunc = [2,2]
                   )

# Plot the sample path 
function plot_tip(x_data, y_data, treshold)
        # Up until the tipping
        lines!(ax1, x_data[1:treshold], y_data[1:treshold], color = (:black,1.00), linewidth = 3.0)
        # After the tipping
        lines!(ax1, x_data[(treshold+1):end], y_data[(treshold+1):end], color = (:black,0.50), linewidth = 3.0)
end

#---------------------------------------------#
#    Subseries, distribution and potential    #
#---------------------------------------------#

function plot_solutions(parameters, subseries, detrended, solution, shifted_potential)
        # Figure for the non-stationary subseries 
        fig1, ax1 = makefig(size = [1300,1800],
                            pad = (20,40,20,20), # Order is: left, right, bottom, top 
                            bg_out = :white,
                            box_position = [1,1:2],
                            limits = ((parameters[1], parameters[end]), (minimum(subseries), maximum(subseries))),
                            toggle_lab = [false, true],
                            lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                            toggle_ticks_lab = [false, true],
                            x_ticks = [parameters[1], parameters[end]],
                            y_ticks = [minimum(subseries), maximum(subseries)],
                            ticks_lab_trunc = [1,1]
                           )
        # Plot the subseries and its trend (qse)
        lines!(ax1, parameters, subseries, color = (:black, 0.25), linewidth = 3.0)
        lines!(ax1, parameters, detrended.trend, color = (:black, 1.0), linewidth = 5.0)

        # Figure for the stationary residuals 
        fig2, ax2 = makefig(fig = fig1,
                            box_position = [2,1:2],
                            limits = ((parameters[1], parameters[end]), (minimum(detrended.residuals), maximum(detrended.residuals))),
                            lab = [L"\mathbf{t}", L"\textbf{residuals}"],
                            lab_pad = [0.0, -40.0],
                            x_ticks = [parameters[1], parameters[end]],
                            y_ticks = [minimum(detrended.residuals), maximum(detrended.residuals)],
                            ticks_lab_trunc = [1,1]
                           )
        # Plot the residuals 
        lines!(ax2, parameters, detrended.residuals, color = (:black, 1.0), linewidth = 3.0)

        # Figure for the distribution 
        fig3, ax3 = makefig(fig = fig1,
                            box_position = [3,1],
                            limits = ((minimum(detrended.residuals), maximum(detrended.residuals)), (0,4)),
                            lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                            x_ticks = [minimum(detrended.residuals), maximum(detrended.residuals)],
                            y_ticks = [0,2,4],
                            ticks_lab_trunc = [1,0]
                           )
        # Compute and plot the histogram of the residuals 
        local bins, pdf = fit_distribution(detrended.residuals, n_bins=convert(Int64, floor(0.005*length(detrended.residuals))))
        domain = LinRange(minimum(bins),maximum(bins),1000)
        barplot!(ax3, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Plot the solution of the NLLS regression to the histogram
        N = get_normalisation_constant(p, solution)
        ρ(x) = N*p(x, solution)
        lines!(ax3, domain, [ρ(x) for x in domain], color = [ρ(x) for x in domain], colormap = [CtpYellow, CtpMauve], linewidth = 8)

        # Figure for the potential 
        fig4, ax4 = makefig(fig = fig1,
                            box_position = [3,2],
                            limits = ((0, 3.5), (0, 3)),
                            lab = [L"\mathbf{x_t}", L"\textbf{V}"],
                            x_ticks = [0, 3.5],
                            y_ticks = [0, 3.0],
                            ticks_lab_trunc = [1,0]
                           )
        # Define the plotting domain for the potential
        domain = LinRange(0,3.5,1000)
        # Loop over the paramer values in the windowed subseries
        for n in 1:length(parameters)
                # Plot only every 100 values of the parameter
                if mod(n, 1000) == 0
                        # Determine transparency value
                        alpha = convert(Float64, n/length(parameters))
                        # Plot the ground truth
                        lines!(ax4, domain, [U(x,parameters[n]) for x in domain], color = (:brown, alpha), linewidth = 1.0)
                end
        end
        # Plot the reconstruction
        lines!(ax4, domain, [shifted_potential(x) for x in domain], color = (:black, 1.0), linewidth = 3.0)
        
        # Adjust whitespace between contiguous plots 
        colgap!(fig1.layout, 40)
        rowgap!(fig1.layout, 50)

        # Return the figure
        return fig1
end

#---------------------------#
#    Solutions and error    #
#---------------------------#

# Figure for the reconstruction error 
fig3, ax3 = makefig(size = [2400, 2000],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    #limits = ((-0.5,8), (0,2)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{||V-U||_2}"],
                    toggle_lab = [false, true],
                    toggle_ticks_lab = [false, true],
                    #x_ticks = [0,4,8],
                    #y_ticks = [0,0.5,1,1.5,2],
                    ticks_lab_trunc = [1,1]
                   )

# Figure for the c1 coefficient 
fig4, ax4 = makefig(fig = fig3,
                    box_position = [1,2],
                    #limits = ((-0.5,8), (0,2)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{c_1}"],
                    toggle_lab = [false, true],
                    toggle_ticks_lab = [false, true],
                    #x_ticks = [0,4,8],
                    #y_ticks = [0,0.5,1,1.5,2],
                    ticks_lab_trunc = [1,1]
                   )

# Figure for the c2 coefficient 
fig5, ax5 = makefig(fig = fig3,
                    box_position = [2,1],
                    #limits = ((-0.5,8), (0,2)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{c_2}"],
                    #x_ticks = [0,4,8],
                    #y_ticks = [0,0.5,1,1.5,2],
                    ticks_lab_trunc = [1,1]
                   )

# Figure for the c3 coefficient 
fig6, ax6 = makefig(fig = fig3,
                    box_position = [2,2],
                    #limits = ((-0.5,8), (0,2)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{c_3}"],
                    #x_ticks = [0,4,8],
                    #y_ticks = [0,0.5,1,1.5,2],
                    ticks_lab_trunc = [1,1]
                   )

function plot_coeffs(timesteps, solutions, error)
        # Define the domain for the x-axis
        domain = LinRange(timesteps[1][1], timesteps[end][end], length(error)) 

        # Extract the coefficients timeseries
        mat = transpose(reduce(hcat, solutions))        # Turns a vector of vectors into a matrix
        c1, c2, c3 = eachcol(mat)                       # Turns the columns of the matrix into vectors
        
        # Plot the error timeseries
        lines!(ax3, domain, error, color = :brown2, linewidth = 3.0) 
        # Plot the c1 timeseries
        lines!(ax4, domain, c1, color = :brown2, linewidth = 3.0) 
        # Plot the c2 timeseries
        lines!(ax5, domain, c2, color = :brown2, linewidth = 3.0) 
        # Plot the c3 timeseries
        lines!(ax6, domain, c3, color = :brown2, linewidth = 3.0) 
end

#------------------------#
#    Random variables    #
#------------------------#

# Figure for the stable equilibrium 
fig7, ax7 = makefig(size = [2400,2000],
                    box_position = [1,1],
                    #limits = ((-0.15,0.15), (0,40)),
                    lab = [L"\mathbf{x_s}", L"\textbf{dist}"],
                    #x_ticks = [-0.15,0,0.15],
                    #y_ticks = [0,20,40],
                    ticks_lab_trunc = [2,2]
                   )

# Figure for the potential value at equilibrium 
fig8, ax8 = makefig(fig = fig7,
                    box_position = [1,2],
                    #limits = ((-0.01,0.001), (0,2000)),
                    lab = [L"\mathbf{V(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    #x_ticks = [-0.01,0],
                    #y_ticks = [0,1000,2000],
                    ticks_lab_trunc = [2,2]
                   )

# Figure for the curvature at equilibrium 
fig9, ax9 = makefig(fig = fig7,
                    box_position = [2,1],
                    #limits = ((0,10), (0,1)),
                    lab = [L"\mathbf{V^{''}(x_s)}", L"\textbf{dist}"],
                    #toggle_lab = [true,false],
                    #x_ticks = [0,10],
                    #y_ticks = [0,0.5,1],
                    ticks_lab_trunc = [2,2]
                   )

# Figure for the exponential of the potential at equilibrium 
fig10, ax10 = makefig(fig = fig7,
                      box_position = [2,2],
                      #limits = ((0,1.05), (0,30)),
                      lab = [L"\textbf{exp}\mathbf{(V(x_s))}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      #x_ticks = [0,1],
                      #y_ticks = [0,10,20,30],
                      ticks_lab_trunc = [2,2]
                     )

function plot_rvs(parameters)
        # Extract the coefficients timeseries
        mat = transpose(reduce(hcat, parameters))       # Turns a vector of vectors into a matrix
        xs, Vxs, Vxxs, eVxs = eachcol(mat)              # Turns the columns of the matrix into vectors

        # Compute the number of bins for the timeseries 
        n_bins = convert(Int64, floor(0.02*length(xs)))

        # Fit and plot a histogram to the stable equilibria 
        bins, pdf = fit_distribution(collect(xs), n_bins=n_bins)
        barplot!(ax7, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the potential values 
        bins, pdf = fit_distribution(Vxs, n_bins=n_bins)
        barplot!(ax8, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the curvature values 
        bins, pdf = fit_distribution(Vxxs, n_bins=n_bins)
        barplot!(ax9, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the exponential of the potential values 
        bins, pdf = fit_distribution(eVxs, n_bins=n_bins)
        barplot!(ax10, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

#-----------------------------#
#    Early-warning signals    #
#-----------------------------#

# Figure for the prefactor 
fig11, ax11 = makefig(size = [2000,1500],
                    box_position = [1,1:2],
                    #limits = ((-0.15,0.15), (0,40)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{C}"],
                    #x_ticks = [-0.15,0,0.15],
                    #y_ticks = [0,20,40],
                    ticks_lab_trunc = [2,4]
                   )

# Figure for the large deviation 
fig12, ax12 = makefig(fig = fig11,
                    box_position = [2,1:2],
                    #limits = ((-0.01,0.001), (0,2000)),
                    lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(\frac{\Delta V}{D})}"],
                    #toggle_lab = [true,false],
                    #x_ticks = [-0.01,0],
                    #y_ticks = [0,1000,2000],
                    ticks_lab_trunc = [2,4]
                   )

# Figure for the escape EWS 
fig13, ax13 = makefig(size = [2000,1500],
                    box_position = [1,1],
                    #limits = ((-0.15,0.15), (0,40)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{C}\textbf{exp}\mathbf{(\frac{\Delta V}{D})}"],
                    #x_ticks = [-0.15,0,0.15],
                    #y_ticks = [0,20,40],
                    ticks_lab_trunc = [2,4]
                   )

function plot_ews(timesteps, escapes)
        # Extract the coefficients timeseries
        mat = transpose(reduce(hcat, escapes))                                      # Turns a vector of vectors into a matrix
        true_prefactor, true_LDP, approx_prefactor, approx_LDP = eachcol(mat)       # Turns the columns of the matrix into vectors

        # Define the domain for the x-axis
        domain = LinRange(timesteps[1][1], timesteps[end][end], length(true_prefactor)) 

        # Plot the prefactor
        lines!(ax11, domain, true_prefactor, color = :brown2, linewidth = 3.0)
        lines!(ax11, domain, approx_prefactor, color = (:black, 0.5), linewidth = 3.0)
        # Plot the large deviation 
        lines!(ax12, domain, true_LDP, color = :brown2, linewidth = 3.0)
        lines!(ax12, domain, approx_LDP, color = (:black, 0.5), linewidth = 3.0)
        # Plot the escape EWS 
        lines!(ax13, domain, true_prefactor.*true_LDP, color = :brown2, linewidth = 3.0)
        lines!(ax13, domain, approx_prefactor.*approx_LDP, color = (:black, 0.5), linewidth = 3.0)
end
