"""
    Figures layout

Generation of the layouts for the figures of the simulations.
"""

#----------------#
#    Figure 1    #
#----------------#

# Figure for the linear reconstructed pdf 
fig1, ax1 = makefig(size = [1600,1600],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-0.35, 0.35), (-0.5, 7)),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    title = L"\mathbf{Nb=%$Nb}",
                    toggle_title = true,
                    x_ticks = [-0.3,0,0.3],
                    y_ticks = [0,7],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the nonlinear reconstructed pdf 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,2],
                    limits = ((-0.35, 0.35), (-0.5, 7)),
                    lab = [L"\mathbf{x_t}", L"\textbf{dist}"],
                    toggle_lab = [true, false],
                    toggle_ticks_lab = [true, false],
                    x_ticks = [-0.3,0,0.3],
                    y_ticks = [0,7],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the linear potential reconstruction 
fig3, ax3 = makefig(fig = fig1,
                    box_position = [2,1],
                    limits = ((-1.5, 1.5), (-1, 1)),
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    x_ticks = [-1.5,0,1.5],
                    y_ticks = [-1,0,1],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the nonlinear potential reconstruction 
fig4, ax4 = makefig(fig = fig1,
                    box_position = [2,2],
                    limits = ((-1.5, 1.5), (-1, 1)),
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    toggle_lab = [true, false],
                    toggle_ticks_lab = [true, false],
                    x_ticks = [-1.5,0,1.5],
                    y_ticks = [-1,0,1],
                    ticks_lab_trunc = [1,0]
                   )

# Adjust whitespace between rows
rowgap!(fig1.layout, 50)
colgap!(fig1.layout, 90)
