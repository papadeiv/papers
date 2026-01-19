"""
    Figures layout
"""

# Figure for the timeseries 
fig, ax = makefig(size = [1200,1200],
                    #pad = (20,40,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    lab = [L"\mathbf{e_{m,t+1}}", L"\mathbf{e_{w,t+1}}"],
                    limits = ((-2, 3.5,), (-2, 3.5)),
                    x_ticks = [-2, 0, 3],
                    y_ticks = [-2, 0, 3],
                    ticks_lab_trunc = [1,1],
                   )
