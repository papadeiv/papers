"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 1200))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Bifurcation diagram
ax1 = Axis(fig[1,1], limits = (-1, 0.025, -1, 1),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [-1; μ_set; 0],
           yticks = [-1,0,1],
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"x",
           ylabelsize = labels,
           yticklabelsize = labels,
          )

# Variance EWS 
ax2 = Axis(fig[2,1], limits = (-1, 0.025, 0, 0.02),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [-1; μ_set; 0],
           yticks = [0,0.02],
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"\text{Var}(x)",
           ylabelsize = labels,
           yticklabelsize = labels,
          )

# Escape EWS 
ax3 = Axis(fig[3,1], limits = (-1, 0.025, 0.25, 1),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [-1; μ_set; 0],
           yticks = [0.25,1],
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\exp(-\Delta V)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )
