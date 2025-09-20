"""
    Figures layout

Generation of the layouts for the figures of the simulations.
"""

CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"
CtpYellow = colorant"rgb(229,200,144)"
CtpWhite = colorant"rgb(198,208,245)"
CtpGray = colorant"rgb(98,104,128)"

#--------------------------------#
#    Timeseries and potential    #
#--------------------------------#

# Figure for the slow-fast timeseries 
fig1, ax1 = makefig(size = [1600,700],
                    pad = (20,50,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1:2],
                    lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                    toggle_lab = [false, true],
                    toggle_ticks_lab = [false, true],
                    ticks_lab_trunc = [3,1]
                   )

# Figure for the detrended residuals 
fig2, ax2 = makefig(fig = fig1,
                    box_position = [2,1:2],
                    lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                    ticks_lab_trunc = [3,1]
                   )

# Figure for the potential 
fig3, ax3 = makefig(fig = fig1,
                    box_position = [1:2,3],
                    limits = ((0.0, 3.5), (-3, 3.5)),
                    lab = [L"\mathbf{x}", L"\mathbf{V}"],
                    x_ticks = [0,3.5],
                    y_ticks = [-3,3.5],
                    ticks_lab_trunc = [1,1]
                   )

# Adjust whitespace between contiguous plots 
colgap!(fig1.layout, 30)
rowgap!(fig1.layout, 40)

#----------------------#
#    NLLS solutions    #
#----------------------#

# Figure for the c1 coefficient
fig4, ax4 = makefig(size = [2100,700],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{c_1}", L"\textbf{dist}"],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the c2 coefficient
fig5, ax5 = makefig(fig = fig4,
                    box_position = [1,2],
                    lab = [L"\mathbf{c_2}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [0,1]
                   )

# Figure for the c3 coefficient
fig6, ax6 = makefig(fig = fig4,
                    box_position = [1,3],
                    lab = [L"\mathbf{c_3}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [0,1]
                   )

#------------------------#
#    Random variables    #
#------------------------#

# Figure for the stable equilibrium 
fig7, ax7 = makefig(size = [2000,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{x_s}", L"\textbf{dist}"],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the unstable equilibrium 
fig8, ax8 = makefig(fig = fig7,
                    box_position = [2,1],
                    lab = [L"\mathbf{x_u}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [0,2]
                   )

# Figure for the potential value at stable equilibrium 
fig9, ax9 = makefig(fig = fig7,
                    box_position = [1,2],
                    lab = [L"\mathbf{V(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the potential value at unstable equilibrium 
fig10, ax10 = makefig(fig = fig7,
                      box_position = [2,2],
                      lab = [L"\mathbf{V(x_u)}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [0,3]
                     )

# Figure for the curvature at stable equilibrium 
fig11, ax11 = makefig(fig = fig7,
                      box_position = [1,3],
                      lab = [L"\mathbf{V^{''}(x_s)}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [0,1]
                     )

# Figure for the curvature at unstable equilibrium 
fig12, ax12 = makefig(fig = fig7,
                      box_position = [2,3],
                      lab = [L"\mathbf{V^{''}(x_u)}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [0,1]
                     )

# Figure for the exponential of the potential at stable equilibrium 
fig13, ax13 = makefig(fig = fig7,
                      box_position = [1,4],
                      lab = [L"\textbf{exp}\mathbf{(V(x_s))}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [2,0]
                     )

# Figure for the exponential of the potential at equilibrium 
fig14, ax14 = makefig(fig = fig7,
                      box_position = [2,4],
                      lab = [L"\textbf{exp}\mathbf{(V(x_u))}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [2,2]
                     )

# Adjust whitespace between contiguous plots 
colgap!(fig4.layout, 50)
rowgap!(fig4.layout, 50)

#-----------------------#
#    Large-deviation    #
#-----------------------#

# Figure for the energy barrier 
fig15, ax15 = makefig(size = [2000,1000],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{\Delta V}", L"\textbf{dist}"],
                    ticks_lab_trunc = [0,2]
                   )

# Figure for the unstable equilibrium 
fig16, ax16 = makefig(fig = fig15,
                    box_position = [1,2],
                    lab = [L"\textbf{exp}\mathbf{(\frac{\Delta V}{D})}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [1,0]
                   )

#=
#------------------#
#    Escape EWS    #
#------------------#

# Figure for the escape ews 
fig11, ax11 = makefig(size = [2000,1000],
                      pad = (30,35,20,20), # Order is: left, right, bottom, top 
                      bg_out = :white,
                      #limits = ((μ0[1],μf[end]), nothing),
                      title = L"\textbf{N_e = %$Ne}", 
                      toggle_title = true,
                      lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(\frac{\Delta U}{D})}"],
                      x_ticks = [μ0[1],μf[end]],
                      #y_ticks = [0,20,40],
                      ticks_lab_trunc = [1,4]
                     )
=#
