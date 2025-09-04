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

# Figure for the timeseries 
fig1, ax1 = makefig(size = [1600,700],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    title = L"\textbf{N_t = %$Nt,\;N_e = %$Ne}",
                    toggle_title = true,
                    limits = ((0, nothing), (-0.3, 0.3)),
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    y_ticks = [-0.3,0.3],
                    ticks_lab_trunc = [0,1]
                   )

# Plot the timeseries 
function plot_ax1(x_data, y_data)
        lines!(ax1, x_data, y_data, color = (:black,0.15), linewidth = 3.0)
end

# Figure for the potential 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,3],
                    limits = ((0.0, 3.5), (0.0, 3.0)),
                    lab = [L"\mathbf{x}", L"\mathbf{V}"],
                    x_ticks = [0,3.5],
                    y_ticks = [0,3],
                    ticks_lab_trunc = [1,0]
                   )
# Plot the ground truth
domain = LinRange(0,3.5,1000)
lines!(ax2, domain, [U(x,μ) for x in domain], color = (:red,1.0), linewidth = 3.0)

# Plot the reconstructed potential 
function plot_ax2(potential)
        lines!(ax2, domain, [potential(x) for x in domain], color = (:black,0.15), linewidth = 1.0)
end

#----------------------#
#    NLLS solutions    #
#----------------------#

# Figure for the c1 coefficient
fig3, ax3 = makefig(size = [2000,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.3,0.3), (0, 12)),
                    lab = [L"\mathbf{c_1}", L"\textbf{dist}"],
                    x_ticks = [-0.3,0.0,0.3],
                    y_ticks = [0,6,12],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the c2 coefficient
fig4, ax4 = makefig(fig = fig3,
                    box_position = [1,2],
                    limits = ((0, 5), (0, 2.5)),
                    lab = [L"\mathbf{c_2}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [0,2.5,5],
                    y_ticks = [0,0.5,1.25,2,2.5],
                    ticks_lab_trunc = [1,1]
                   )

# Figure for the c3 coefficient
fig5, ax5 = makefig(fig = fig3,
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
        n_bins = convert(Int64, floor(0.1*Ne))

        # Fit and plot a histogram to the c1 coefficients
        bins, pdf = fit_distribution(coefficients[:,1], n_bins=n_bins)
        barplot!(ax3, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the c2 coefficients
        bins, pdf = fit_distribution(coefficients[:,2], n_bins=n_bins)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the c3 coefficients
        bins, pdf = fit_distribution(coefficients[:,3], n_bins=n_bins)
        barplot!(ax5, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

#-------------#
#    Error    #
#-------------#

# Figure for the reconstruction error 
Fig, Ax = makefig(fig = fig3,
                  box_position = [1,4],
                  limits = ((-0.5,8), (0,2)),
                  lab = [L"\mathbf{||V-U||_2}", L"\textbf{dist}"],
                  toggle_lab = [true,false],
                  x_ticks = [0,4,8],
                  y_ticks = [0,0.5,1,1.5,2],
                  ticks_lab_trunc = [0,1]
                 )

function plot_error(error)
        # Compute the number of bins for the ensemble
        n_bins = convert(Int64, floor(0.02*Ne))

        # Fit and plot a histogram to the ensemble reconstruction error 
        bins, pdf = fit_distribution(error, n_bins=n_bins)
        barplot!(Ax, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

#------------------------#
#    Random variables    #
#------------------------#

# Figure for the stable equilibrium 
fig6, ax6 = makefig(fig = fig3,
                    box_position = [2,1],
                    limits = ((-0.15,0.15), (0,40)),
                    lab = [L"\mathbf{x_s}", L"\textbf{dist}"],
                    x_ticks = [-0.15,0,0.15],
                    y_ticks = [0,20,40],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the potential value at equilibrium 
fig7, ax7 = makefig(fig = fig3,
                    box_position = [2,2],
                    limits = ((-0.01,0.001), (0,2000)),
                    lab = [L"\mathbf{V(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [-0.01,0],
                    y_ticks = [0,1000,2000],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the curvature at equilibrium 
fig8, ax8 = makefig(fig = fig3,
                    box_position = [2,3],
                    limits = ((0,10), (0,1)),
                    lab = [L"\mathbf{V^{''}(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [0,10],
                    y_ticks = [0,0.5,1],
                    ticks_lab_trunc = [0,1]
                   )

# Figure for the exponential of the potential at equilibrium 
fig9, ax9 = makefig(fig = fig3,
                    box_position = [2,4],
                    limits = ((0,1.05), (0,30)),
                    lab = [L"\textbf{exp}\mathbf{(V^{''}(x_s))}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    x_ticks = [0,1],
                    y_ticks = [0,10,20,30],
                    ticks_lab_trunc = [0,0]
                   )

function plot_rvs(variables)
        # Compute the number of bins for the ensemble
        n_bins = convert(Int64, floor(0.1*Ne))

        # Fit and plot a histogram to the stable equilibria 
        bins, pdf = fit_distribution(variables[:,1], n_bins=n_bins)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the potential values 
        bins, pdf = fit_distribution(variables[:,2], n_bins=n_bins)
        barplot!(ax7, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the curvature values 
        bins, pdf = fit_distribution(variables[:,3], n_bins=n_bins)
        barplot!(ax8, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)

        # Fit and plot a histogram to the exponential of the potential values 
        bins, pdf = fit_distribution(variables[:,4], n_bins=n_bins)
        barplot!(ax9, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
end

# Adjust whitespace between contiguous plots 
colgap!(fig3.layout, 30)
rowgap!(fig3.layout, 30)

# Insert title
Label(fig3[begin-1, 1:4],
      L"\textbf{N_t = %$Nt,\; N_e = %$Ne}",
      fontsize = 50,
      padding = (0,0,0,0),
)
