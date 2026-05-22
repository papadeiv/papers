"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1800, 600), 
             #figure_padding = (60,60,30,30), 
             #backgroundcolor = :white
            )

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Linear coefficient
ax1 = Axis(fig[1,1], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = c,
           #yticks = [0,9],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\alpha",
           ylabel = L"\text{distribution}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
           titlevisible = true,
           titlesize = labels,
          )

# Mean 
ax2 = Axis(fig[1,2], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = c,
           #yticks = [0,9],
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

# Standard deviation 
ax3 = Axis(fig[1,3], limits = (nothing, nothing, 0, nothing),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = false,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           #xticks = c,
           #yticks = [0,9],
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

# Array of axes
axes = [ax1, ax2, ax3]

# Adjust whitespace between columns and rows
#colgap!(fig.layout, 120)
#rowgap!(fig.layout, 10)
