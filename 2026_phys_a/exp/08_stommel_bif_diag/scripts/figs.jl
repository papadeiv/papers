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

# Bifurcation diagram in x (ψ)
ax1 = Axis(fig[1,1], limits = (0, 2, -1, 1.5),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           yticks = [-1,1.5],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\psi",
           xlabelsize = labels,
           ylabelsize = labels,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Bifurcation diagram in y (S)
ax2 = Axis(fig[1,2], limits = (0, 2, -1, 3.5),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           yticks = [-1,3.5],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"S",
           xlabelsize = labels,
           ylabelsize = labels,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
