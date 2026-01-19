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

#----------------#
#    Figure 2    #
#----------------#

# Figure for the c1 coefficient
fig3, ax3 = makefig(size = [1500,500],
                    pad = (30,35,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{c_1}", L"\textbf{dist}"],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for the c2 coefficient
fig4, ax4 = makefig(fig = fig3,
                    box_position = [1,2],
                    lab = [L"\mathbf{c_2}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for the c3 coefficient
fig5, ax5 = makefig(fig = fig3,
                    box_position = [1,3],
                    lab = [L"\mathbf{c_3}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [1,0]
                   )

rowgap!(fig3.layout, 50)

#----------------#
#    Figure 3    #
#----------------#

# Figure for xs 
fig6, ax6 = makefig(size = [3000,1500],
                    pad = (30,50,20,20), # Order is: left, right, bottom, top 
                    bg_out = :white,
                    box_position = [1,1],
                    lab = [L"\mathbf{x_s}", L"\textbf{dist}"],
                    ticks_lab_trunc = [3,0]
                   )

# Figure for V(xs) 
fig7, ax7 = makefig(fig = fig6,
                    box_position = [1,2],
                    lab = [L"\mathbf{V(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [2,0]
                   )

# Figure for V''(xs) 
fig8, ax8 = makefig(fig = fig6,
                    box_position = [1,3],
                    lab = [L"\mathbf{V^{''}(x_s)}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [1,0]
                   )

# Figure for e^(V(xs))
fig9, ax9 = makefig(fig = fig6,
                    box_position = [1,4],
                    lab = [L"\textbf{exp}\mathbf{(V(x_s))}", L"\textbf{dist}"],
                    toggle_lab = [true,false],
                    ticks_lab_trunc = [3,0]
                   )

# Figure for xu 
fig10, ax10 = makefig(fig = fig6,
                      box_position = [2,1],
                      lab = [L"\mathbf{x_u}", L"\textbf{dist}"],
                      ticks_lab_trunc = [1,1]
                     )

# Figure for V(xu) 
fig11, ax11 = makefig(fig = fig6,
                      box_position = [2,2],
                      lab = [L"\mathbf{V(x_u)}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [1,1]
                     )

# Figure for V''(xu) 
fig12, ax12 = makefig(fig = fig6,
                      box_position = [2,3],
                      lab = [L"\mathbf{V^{''}(x_u)}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [1,0]
                     )

# Figure for e^(V(xu))
fig13, ax13 = makefig(fig = fig6,
                      box_position = [2,4],
                      lab = [L"\textbf{exp}\mathbf{(V(x_u))}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [1,1]
                     )

rowgap!(fig6.layout, 50)

#----------------#
#    Figure 4    #
#----------------#

# Figure for ΔV
fig14, ax14 = makefig(size = [1500,1500],
                      pad = (30,35,20,20), # Order is: left, right, bottom, top 
                      bg_out = :white,
                      box_position = [1,1],
                      lab = [L"\mathbf{ΔV}", L"\textbf{dist}"],
                      ticks_lab_trunc = [1,0]
                     )

# Figure for sqrt(|V''(xu)|V''(xs))
fig15, ax15 = makefig(fig = fig14,
                      box_position = [1,2],
                      lab = [L"\mathbf{\sqrt{|V^{''}(x_u)|V^{''}(x_s)}}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [1,0]
                     )

# Figure for exp(-ΔV)
fig16, ax16 = makefig(fig = fig14,
                      box_position = [2,1],
                      lab = [L"\textbf{exp}\mathbf{(-ΔV)}", L"\textbf{dist}"],
                      ticks_lab_trunc = [2,0]
                     )

# Figure for the numerical error
fig17, ax17 = makefig(fig = fig14,
                      box_position = [2,2],
                      lab = [L"\mathbf{||V-V_{*}||_2}", L"\textbf{dist}"],
                      toggle_lab = [true,false],
                      ticks_lab_trunc = [2,0]
                     )

rowgap!(fig14.layout, 50)
colgap!(fig14.layout, 50)
