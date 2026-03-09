"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [2000,1000],
                    pad = (20,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-2.5,2.5),(-0.01,0.6)),
                    lab = [L"\mathbf{x}", L"\textbf{distribution}"],
                    x_ticks = [-2,0,2],
                    y_ticks = [0,0.5],
                    ticks_lab_trunc = [0,1]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,2],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    limits = ((-1.5, 1.5), (-0.05, 1.5)),
                    x_ticks = [-1, 0, 1],
                    y_ticks = [0, 1],
                    ticks_lab_trunc = [0,0]
                   )
