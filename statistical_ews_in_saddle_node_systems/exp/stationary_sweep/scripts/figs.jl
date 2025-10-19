"""
    Figures layout

Generation of the layouts for the figures of the simulations.
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

#----------------#
#    Figure 2    #
#----------------#

# Figure for the large deviation principle (LDP)
fig3, ax3 = makefig(size = [1500,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-1.05,-0.05),(-0.05,1.05)),
                    lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(-ΔV)}"],
                    x_ticks = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1],
                    y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [1,2]
                   )

#----------------#
#    Figure 3    #
#----------------#

# Figure for the prefactor (C)
fig4, ax4 = makefig(size = [1500,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-1.05,-0.05),(nothing,nothing)),
                    lab = [L"\mathbf{\mu}", L"\mathbf{\sqrt{|V^{''}(b)|V^{''}(a)}}"],
                    x_ticks = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1],
                    #y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [1,2]
                   )

#----------------#
#    Figure 4    #
#----------------#

# Figure for the early-warning signal (EWS)
fig5, ax5 = makefig(size = [1500,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    limits = ((-1.05,-0.05),(nothing,nothing)),
                    lab = [L"\mathbf{\mu}", L"\textbf{ews}"],
                    x_ticks = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1],
                    #y_ticks = [0,0.25,0.5,0.75,1],
                    ticks_lab_trunc = [1,2]
                   )
