"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#----------------#
#    Figure 1    #
#----------------#

# Figure for the timeseries 
fig1, ax1 = makefig(size = [1600,600],
                    pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    title = L"\textbf{N_t = %$Nt,\;N_e = %$Ne}",
                    toggle_title = false,
                    lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                    ticks_lab_trunc = [0,1]
                   )

# Figure for the potential 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [1,3],
                    limits = ((-1.5, 1.5), (-1, 1)),
                    lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                    x_ticks = [-1.5,0,1.5],
                    y_ticks = [-1,0,1],
                    ticks_lab_trunc = [1,0]
                   )
