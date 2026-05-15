"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

fig = Figure(; size = (900,1200))
ax1 = Axis(fig[1:2,1:2], limits = (c0, 3, nothing, nothing))
ax2 = Axis(fig[3:4,1:2], limits = (c0, 3, nothing, nothing))
ax3L = Axis(fig[5,1], title = L"\textbf{variance}")
ax3R = Axis(fig[5,2], title = L"\textbf{skewness}")
ax4L = Axis(fig[6,1], title = L"\textbf{autocorr}")
ax4R = Axis(fig[6,2], title = L"\textbf{exit rate}")
