"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

fig = Figure(; size = (900,1200))
ax1 = Axis(fig[1:2,1:2], limits = (c0, cf, nothing, nothing))
ax2 = Axis(fig[3:4,1:2], limits = (c0, cf, nothing, nothing))
ax3L = Axis(fig[5,1])
ax3R = Axis(fig[5,2])
ax4L = Axis(fig[6,1])
ax4R = Axis(fig[6,2])
