"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1300, 650), figure_padding = (30,30,30,30))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Variance EWS 
ax1 = Axis(fig[1,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=3))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\text{Var}(x)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Escape EWS 
ax2 = Axis(fig[1,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\exp(-\Delta V)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
