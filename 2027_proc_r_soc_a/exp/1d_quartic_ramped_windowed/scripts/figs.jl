"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

fig1, ax1 = makefig(size = [1600,1600],
                    pad = (20,60,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    toggle_lab = [false,true],
                    toggle_ticks_lab = [false,true],
                    ticks_lab_trunc = [1,2]
                   )

fig2, ax2 = makefig(fig = fig1,
                    box_position = [2,1],
                    lab = [L"\mathbf{t}", L"\textbf{exp}\mathbf{(-\Delta V)}"],
                    toggle_lab = [false,true],
                    toggle_ticks_lab = [false,true],
                    ticks_lab_trunc = [1,2]
                   )

fig3, ax3 = makefig(fig = fig1,
                    box_position = [3,1],
                    lab = [L"\mathbf{\mu}", L"\mathbf{||V - V_{*}||_2}"],
                    ticks_lab_trunc = [1,2]
                   )

rowgap!(fig1.layout, 50)

#----------------#
#    Figure 2    #
#----------------#

fig4, ax4 = makefig(size = [2000,1000],
                    pad = (20,60,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1:2,1:2],
                    lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                    toggle_ticks = [false,true],
                    ticks_lab_trunc = [1,1]
                   )

fig5, ax5 = makefig(fig = fig4,
                    box_position = [1,3],
                    lab = [L"\mathbf{x}", L"\mathbf{p(x;\,\mu)}"],
                    y_ticks = [0, 16],
                    ticks_lab_trunc = [1,1]
                   )

fig6, ax6 = makefig(fig = fig4,
                    box_position = [2,3],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x;\,\mu)}"],
                    limits = ((-1.5, 1.5), (-2, 2)),
                    x_ticks = [-1.5,0,1.5],
                    y_ticks = [-2,2],
                    ticks_lab_trunc = [1,1]
                   )
