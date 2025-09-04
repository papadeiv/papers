"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"
CtpYellow = colorant"rgb(229,200,144)"

#--------------------------------#
#    Timeseries and potential    #
#--------------------------------#

# Figure for the slow-fast timeseries 
fig1, ax1 = makefig(size = [1600,700],
                    pad = (20,50,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    limits = ((μ0, μf), (2.3, 3.1)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                    toggle_lab = [false, true],
                    toggle_ticks_lab = [false, true],
                    x_ticks = [μ0, μf],
                    y_ticks = [2.3, 3.1],
                    ticks_lab_trunc = [3,1]
                   )

# Figure for the detrended residuals 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [2,1:2],
                    limits = ((μ0, μf), (-0.4, 0.4)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                    x_ticks = [μ0, μf],
                    y_ticks = [-0.4,0.4],
                    ticks_lab_trunc = [1,1]
                   )

# Plot the timeseries 
function plot_solutions(x_data, full, residual)
        lines!(ax1, x_data, full, color = (:black,0.15), linewidth = 3.0)
        lines!(ax2, x_data, residual, color = (:black,0.15), linewidth = 3.0)
end

# Plot the drift of the QSE 
function plot_drift(x_data, trend)
        lines!(ax1, x_data, trend, color = (CtpYellow,1.00), linewidth = 5.0)
end

# Figure for the potential 
fig3, ax3 = makefig(fig = fig1,
                    box_position = [1:2,3],
                    limits = ((0.0, 3.5), (-1, 1)),
                    lab = [L"\mathbf{x}", L"\mathbf{V}"],
                    x_ticks = [0,3.5],
                    y_ticks = [-1,1],
                    ticks_lab_trunc = [1,0]
                   )
# Define the domain for plotting the potential 
domain = LinRange(0,3.5,1000)

# Plot the true scalar potential 
function plot_ground_truth(parameters)
        # Loop over the paramer values in the windowed subseries
        for n in 1:length(parameters)
                # Plot only every 100 values of the parameter
                if mod(n, 1000) == 0
                        # Determine transparency value
                        alpha = convert(Float64, n/length(parameters))
                        # Plot the ground truth
                        lines!(ax3, domain, [U(x,parameters[n]) for x in domain], color = (:red, alpha), linewidth = 1.0)
                end
        end
end

# Plot the reconstructed potential 
function plot_reconstruction(potential)
        lines!(ax3, domain, [potential(x) for x in domain], color = (:black,0.15), linewidth = 1.0)
end

# Adjust whitespace between contiguous plots 
colgap!(fig1.layout, 30)
rowgap!(fig1.layout, 40)

#----------------------#
#    NLLS solutions    #
#----------------------#

# Figure for the c1 coefficient
fig4, ax4 = makefig(size = [2000,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.5,0.5), (0,6)),
                    lab = [L"\mathbf{c_1}", L"\textbf{dist}"],
                    x_ticks = [-0.5,0.5],
                    y_ticks = [0,3,6],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the c2 coefficient
fig5, ax5 = makefig(fig = fig4,
                    box_position = [1,2],
                    limits = ((0,5), (0,2)),
                    lab = [L"\mathbf{c_2}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [0,5],
                    y_ticks = [0,0.5,1.5,2],
                    ticks_lab_trunc = [0,1]
                   )

# Figure for the c3 coefficient
fig6, ax6 = makefig(fig = fig4,
                    box_position = [1,3],
                    limits = ((-20,20), (0,0.3)),
                    lab = [L"\mathbf{c_3}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [-20,-5,5,20],
                    y_ticks = [0,0.15,0.3],
                    ticks_lab_trunc = [0,2]
                   )

function plot_coeffs(coefficients)
        # Compute the number of bins for the ensemble
        n_bins = convert(Int64, floor(0.10*Ne))

        # Fit and plot a histogram to the c1 coefficients
        bins, pdf = fit_distribution(coefficients[:,1], n_bins=n_bins)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the c2 coefficients
        bins, pdf = fit_distribution(coefficients[:,2], n_bins=n_bins)
        barplot!(ax5, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the c3 coefficients
        bins, pdf = fit_distribution(coefficients[:,3], n_bins=n_bins)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

#------------------------#
#    Random variables    #
#------------------------#

# Figure for the stable equilibrium 
fig7, ax7 = makefig(fig = fig4,
                    box_position = [2,1],
                    limits = ((-0.15,0.15), (0,40)),
                    lab = [L"\mathbf{x_s}", L"\textbf{dist}"],
                    x_ticks = [-0.15,0,0.15],
                    y_ticks = [0,40],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the potential value at equilibrium 
fig8, ax8 = makefig(fig = fig4,
                    box_position = [2,2],
                    limits = ((-0.01,0.001), (0,1000)),
                    lab = [L"\mathbf{V(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [-0.01,0],
                    y_ticks = [0,500,1000],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the curvature at equilibrium 
fig9, ax9 = makefig(fig = fig4,
                    box_position = [2,3],
                    limits = ((0,10), (0,1)),
                    lab = [L"\mathbf{V^{''}(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [0,10],
                    y_ticks = [0,0.5,1],
                    ticks_lab_trunc = [0,1]
                   )

# Figure for the exponential of the potential at equilibrium 
fig10, ax10 = makefig(fig = fig4,
                      box_position = [2,4],
                      limits = ((0,1.05), (0,40)),
                      lab = [L"\textbf{exp}\mathbf{(V(x_s))}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      x_ticks = [0,1],
                      y_ticks = [0,10,20,30,40],
                      ticks_lab_trunc = [0,0]
                     )

function plot_rvs(variables)
        # Compute the number of bins for the ensemble
        n_bins = convert(Int64, floor(0.10*Ne))

        # Fit and plot a histogram to the stable equilibria 
        bins, pdf = fit_distribution(variables[:,1], n_bins=n_bins)
        barplot!(ax7, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the potential values 
        bins, pdf = fit_distribution(variables[:,2], n_bins=n_bins)
        barplot!(ax8, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the curvature values 
        bins, pdf = fit_distribution(variables[:,3], n_bins=n_bins)
        barplot!(ax9, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the exponential of the potential values 
        bins, pdf = fit_distribution(variables[:,4], n_bins=n_bins)
        barplot!(ax10, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

# Adjust whitespace between contiguous plots 
colgap!(fig4.layout, 40)
rowgap!(fig4.layout, 30)

# Print the title with info
function print_info(Nt)
        Label(fig1[begin-1, 1:2],
              L"\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
        Label(fig4[begin-1, 1:4],
              L"\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
end
