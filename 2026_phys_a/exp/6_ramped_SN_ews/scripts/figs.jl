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
ax1 = Axis(fig[1,1:2], limits = (-1, 0.025, -1.1, 1.1),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [μ0, μf],
           yticks = [-1,0,1],
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"x",
           ylabelsize = labels,
           yticklabelsize = labels,
          )

# Residuals 
ax2 = Axis(fig[2,1:2], limits = (-1, 0.025, -0.5, 0.5),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = [μ0, μf],
           yticks = [-0.5,0,0.5],
           ytickformat = values -> ["$(trunc(value, digits=2))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"x",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Variance EWS 
ax3_L = Axis(fig[3,1],
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
ax3_R = Axis(fig[3,2],
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
             ytickformat = values -> ["$(trunc(value, digits=0))" for value in values],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"\mu",
             ylabel = L"\exp(-\Delta V)",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )
