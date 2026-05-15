"""
    Figures layout

Generation of the layouts for the figures of the simulations.
"""

#----------------#
#    Figure 1    #
#----------------#

# Figure for the OUP histogram 
fig1, ax1 = makefig(size = [1600,1000],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1:3,1],
                    limits = ((-0.4, 0.4), (-0.5, 7)),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    toggle_lab = [false, true],
                    title = L"\textbf{OUP } \mathbf{\theta=%$θ}",
                    toggle_title = true,
                    x_ticks = [-0.3,0,0.3],
                    y_ticks = [0,7],
                    toggle_ticks_lab = [false, true],
                    ticks_lab_trunc = [1,0]
                   )


# Figure for the Gaussian histogram 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [1:3,2],
                    limits = ((-0.4, 0.4), (-0.5, 7)),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    toggle_lab = [false, false],
                    title = L"\textbf{Gaussian } \mathbf{\sigma=%$(round(sqrt((σ^2)/(2*θ)), digits=3))}",
                    toggle_title = true,
                    toggle_ticks_lab = [false, false],
                    x_ticks = [-0.3,0,0.3],
                    y_ticks = [0,7],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the OUP samples 
fig3, ax3 = makefig(fig = fig1,
                    box_position = [4,1],
                    limits = ((-0.4, 0.4), nothing),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    toggle_lab = [true, false],
                    toggle_ticks = [true, false],
                    toggle_ticks_lab = [true, false],
                    x_ticks = [-0.3,0,0.3],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the Gaussian drawn samples 
fig4, ax4 = makefig(fig = fig1,
                    box_position = [4,2],
                    limits = ((-0.4, 0.4), nothing),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    toggle_lab = [true, false],
                    toggle_ticks = [true, false],
                    toggle_ticks_lab = [true, false],
                    x_ticks = [-0.3,0,0.3],
                    ticks_lab_trunc = [1,0]
                   )

# Adjust whitespace between rows
rowgap!(fig1.layout, 50)
colgap!(fig1.layout, 50)

#----------------#
#    Figure 2    #
#----------------#

# Figure for the OUP histogram 
fig5, ax5 = makefig(size = [1600,2400],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.4, 0.4), nothing),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    title = L"\textbf{OUP } \mathbf{\theta=%$θ}",
                    toggle_title = true,
                    x_ticks = [-0.3,0,0.3],
                    #y_ticks = [0,7],
                    ticks_lab_trunc = [1,2]
                   )


# Figure for the Gaussian histogram 
fig6, ax6 = makefig(fig = fig5,
                    box_position = [2,1],
                    limits = ((-0.4, 0.4), nothing),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    title = L"\textbf{Gaussian } \mathbf{\sigma=%$(round(sqrt((σ^2)/(2*θ)), digits=3))}",
                    toggle_title = true,
                    x_ticks = [-0.3,0,0.3],
                    #y_ticks = [0,7],
                    ticks_lab_trunc = [1,2]
                   )


