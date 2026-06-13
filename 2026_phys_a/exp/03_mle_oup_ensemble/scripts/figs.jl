"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1696, 848))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Linear coefficient (α = 10)
ax1 = Axis(fig[1,1], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\alpha",
           ylabel = L"\text{distribution (deep)}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Mean (α = 10) 
ax2 = Axis(fig[1,2], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"a",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Standard deviation (α = 10)
ax3 = Axis(fig[1,3], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\sqrt{2D}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Linear coefficient (α = 1)
ax4 = Axis(fig[2,1], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\alpha",
           ylabel = L"\text{distribution (shallow)}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Mean (α = 1) 
ax5 = Axis(fig[2,2], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"a",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Standard deviation (α = 1)
ax6 = Axis(fig[2,3], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\sqrt{2D}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Potential (α = 10)
ax_p1 = Axis(fig[1,4], limits = (-0.5, 0.5, -0.01, 0.5),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [-0.5, 0.5],
           yticks = [0, 0.5],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"x",
           ylabel = L"\text{potential (deep)}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Potential (α = 1)
ax_p2 = Axis(fig[2,4], limits = (-0.5, 0.5, -0.01, 0.5),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [-0.5, 0.5],
           yticks = [0, 0.5],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"x",
           ylabel = L"\text{potential (shallow)}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Arrays of axes
axes = [ax1, ax2, ax3, ax4, ax5, ax6]
ax_p = [ax_p1, ax_p2]
