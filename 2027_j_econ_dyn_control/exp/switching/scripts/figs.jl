"""
    Figures layout
"""

#------------------#
#    Time-paths    #
#------------------#

fig, ax1 = makefig(size = [2100,900],
                   #pad = (20,40,20,20), # Order is: left, right, bottom, top 
                   bg_out = :white,
                   #bg_in = :transparent,
                   box_position = [1,1],
                   lab = [L"\mathbf{t}", L"\mathbf{g_{t}}"],
                   toggle_lab = [false, true],
                   limits = ((0,N), (nothing, nothing)),
                   x_ticks = [0,N],
                   y_ticks = [0,2.4],
                   toggle_ticks_lab = [false, true],
                   ticks_lab_trunc = [0,1],
                  )

fig, ax2 = makefig(fig = fig,
                   box_position = [1,2],
                   lab = [L"\mathbf{t}", L"\mathbf{A_{t}}"],
                   toggle_lab = [false, true],
                   limits = ((0,N), (0, 25)),
                   x_ticks = [0,N],
                   y_ticks = [0,25],
                   toggle_ticks_lab = [false, true],
                   ticks_lab_trunc = [0,0],
                  )

fig, ax3 = makefig(fig = fig,
                   box_position = [2,1],
                   lab = [L"\mathbf{t}", L"\mathbf{e_{t}}"],
                   limits = ((0,N), (nothing, nothing)),
                   x_ticks = [0,N],
                   y_ticks = [0,0.07],
                   ticks_lab_trunc = [0,2],
                  )

fig, ax4 = makefig(fig = fig,
                   box_position = [2,2],
                   lab = [L"\mathbf{t}", L"\mathbf{L_{t}}"],
                   limits = ((0,N), (nothing, nothing)),
                   x_ticks = [0,N],
                   y_ticks = [0,16],
                   ticks_lab_trunc = [0,0],
                  )
