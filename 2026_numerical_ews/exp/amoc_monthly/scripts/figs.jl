"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [1400,700],
                    pad = (20,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    limits = ((0,2200),(-7,25)),
                    lab = [L"\textbf{time (millenia)}", L"\textbf{observable (Sv)}"],
                    toggle_lab = [false,true],
                    x_ticks = [0,1000,2000],
                    y_ticks = [-7,25],
                    toggle_ticks_lab = [false,true],
                    ticks_lab_trunc = [0,0]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [1:2,3:4],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    limits = ((-5, 5), (-5, 5)),
                    x_ticks = [-5, 5],
                    y_ticks = [-5, 5],
                    ticks_lab_trunc = [0,0]
                   )

fig3, ax3 = makefig(fig = fig1,
                    box_position = [2,1:2],
                    limits = ((0,2200),(-0.05,1.05)),
                    lab = [L"\textbf{time (millenia)}", L"\textbf{exp}\mathbf{(-ΔV)}"],
                    x_ticks = [0,1000,2000],
                    y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [0,2]
                   )
