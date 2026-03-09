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
                    limits = ((-1,1),(-0.1,3)),
                    lab = [L"\mathbf{x}", L"\textbf{distribution}"],
                    x_ticks = [-1,0,1],
                    y_ticks = [0,3],
                    ticks_lab_trunc = [0,0]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,2],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    limits = ((-1, 1), (0, 1)),
                    x_ticks = [-1, 0, 1],
                    y_ticks = [0, 1],
                    ticks_lab_trunc = [0,0]
                   )
