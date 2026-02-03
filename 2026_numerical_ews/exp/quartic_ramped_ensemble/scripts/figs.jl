"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

fig1, ax1 = makefig(size = [900,900],
                    pad = (20,60,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{x}", L"\mathbf{V(x;\,\mu)}"],
                    limits = ((-1.5, 1.5), (-2, 2)),
                    x_ticks = [-1.5, 0, 1.5],
                    y_ticks = [-2, 0, 2],
                    ticks_lab_trunc = [1,1]
                   )
