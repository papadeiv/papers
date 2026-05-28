"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 600))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

ax1 = Axis(fig[1,1], limits = (-2, 2, -1, 1),
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = true,
             xticklabelsvisible = true,
             yticklabelsvisible = true,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [-2,0,2],
             yticks = [-1,0,1],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

ax2 = Axis(fig[1,2], limits = (nothing, nothing, nothing, nothing),
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
             xlabel = L"\text{cond. number}",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )
