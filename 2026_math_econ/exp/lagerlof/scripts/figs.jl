"""
    Figures layout
"""

#------------------#
#    Time-paths    #
#------------------#

fig, ax1 = makefig(size = [2400,1600],
                   #pad = (20,40,20,20), # Order is: left, right, bottom, top 
                   bg_out = :transparent,
                   bg_in = :transparent,
                   box_position = [1,1:3],
                   lab = [L"\mathbf{t}", L"\mathbf{g_{t}}"],
                   toggle_lab = [false, true],
                   x_ticks = [0, 50, N],
                   toggle_ticks_lab = [false, true],
                   ticks_lab_trunc = [0,1],
                  )

fig, ax2 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [1,4:6],
                   lab = [L"\mathbf{t}", L"\mathbf{A_{t}}"],
                   toggle_lab = [false, true],
                   limits = ((0,300), (-0.05*16.5, 1.05*16.5)),
                   x_ticks = [0, 50, N],
                   y_ticks = [0, 16],
                   toggle_ticks_lab = [false, true],
                   ticks_lab_trunc = [0,0],
                  )

fig, ax3 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [2,1:3],
                   lab = [L"\mathbf{t}", L"\mathbf{e_{t}}"],
                   x_ticks = [0, 50, N],
                   ticks_lab_trunc = [0,2],
                  )

fig, ax4 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [2,4:6],
                   lab = [L"\mathbf{t}", L"\mathbf{L_{t}}"],
                   x_ticks = [0, 50, N],
                   ticks_lab_trunc = [0,0],
                  )

fig, ax5 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [3,1:2],
                   lab = [L"\mathbf{t}", L"\mathbf{g_{t}}"],
                   x_ticks = [0, 50],
                   ticks_lab_trunc = [0,1],
                  )

fig, ax6 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [3,3:4],
                   lab = [L"\mathbf{t}", L"\mathbf{e_{t}}"],
                   x_ticks = [0, 50],
                   ticks_lab_trunc = [0,2],
                  )

fig, ax7 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [3,5:6],
                   lab = [L"\mathbf{t}", L"\mathbf{L_{t}}"],
                   x_ticks = [0, 50],
                   ticks_lab_trunc = [0,0],
                  )

rowgap!(fig.layout, 2, Relative(0.05))

#-------------------#
#    Observables    #
#-------------------#

fig, ax8 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [4,1:2],
                   lab = [L"\mathbf{t}", L"\mathbf{h_{t}}"],
                   limits = ((0,300), (0, 1)),
                   x_ticks = [0, 50, N],
                   y_ticks = [0, 1],
                   ticks_lab_trunc = [0,0],
                  )

fig, ax9 = makefig(fig = fig,
                   bg_in = :transparent,
                   box_position = [4,3:4],
                   lab = [L"\mathbf{t}", L"\mathbf{x_{t}}"],
                   limits = ((0,300), (0, 8)),
                   x_ticks = [0, 50, N],
                   y_ticks = [0, 8],
                   ticks_lab_trunc = [0,0],
                  )

fig, ax10 = makefig(fig = fig,
                    bg_in = :transparent,
                    box_position = [4,5:6],
                    lab = [L"\mathbf{t}", L"\mathbf{z_{t}}"],
                    limits = ((0,300), (0, 5)),
                    x_ticks = [0, 50, N],
                    y_ticks = [0, 5],
                    ticks_lab_trunc = [0,0],
                   )

rowgap!(fig.layout, 3, Relative(0.05))
