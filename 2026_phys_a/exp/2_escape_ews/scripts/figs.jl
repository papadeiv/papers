"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 1200), 
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

# Bifurcation diagram
ax1 = Axis(fig[1:2,1:5], limits = (c0, 3, 0, 9),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = true,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = c,
           yticks = [0,9],
           xtickwidth = border,
           ytickwidth = border,
           ylabel = L"x",
           ylabelsize = labels,
           yticklabelsize = labels,
          )

# Potential at μ = 1.2 
ax2_1 = Axis(fig[3,1], limits = (0, 9, lower, upper),
             title = L"\mu=%$(c[1])",
             titlesize = labels,
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = true,
             xticklabelsvisible = true,
             yticklabelsvisible = true,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [0,9],
             yticks = [lower,upper],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             ylabel = L"V(x\,;\;\mu)",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

# Potential at μ = 1.6 
ax2_2 = Axis(fig[3,2], limits = (0, 9, lower, upper),
             title = L"\mu=%$(c[2])",
             titlesize = labels,
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = false,
             xticklabelsvisible = true,
             yticklabelsvisible = false,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [0,9],
             yticks = [lower,upper],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             xlabelsize = labels,
             xticklabelsize = labels,
            )

# Potential at μ = 2.0 
ax2_3 = Axis(fig[3,3], limits = (0, 9, lower, upper),
             title = L"\mu=%$(c[3])",
             titlesize = labels,
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = false,
             xticklabelsvisible = true,
             yticklabelsvisible = false,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [0,9],
             yticks = [lower,upper],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             xlabelsize = labels,
             xticklabelsize = labels,            
            )

# Potential at μ = 2.4 
ax2_4 = Axis(fig[3,4], limits = (0, 9, lower, upper),
             title = L"\mu=%$(c[4])",
             titlesize = labels,
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = false,
             xticklabelsvisible = true,
             yticklabelsvisible = false,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [0,9],
             yticks = [lower,upper],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             xlabelsize = labels,
             xticklabelsize = labels,            
            )

# Potential at μ = 2.8 
ax2_5 = Axis(fig[3,5], limits = (0, 9, lower, upper),
             title = L"\mu=%$(c[5])",
             titlesize = labels,
             spinewidth = border,
             xgridvisible = false,
             ygridvisible = false,
             xlabelvisible = true,
             ylabelvisible = false,
             xticklabelsvisible = true,
             yticklabelsvisible = false,
             xtickalign = 1,
             ytickalign = 1,
             xticks = [0,9],
             yticks = [lower,upper],
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"x",
             xlabelsize = labels,
             xticklabelsize = labels,            
            )

# Escape EWS 
ax3 = Axis(fig[4:5,1:5], limits = (c0,3,0,1),
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xticks = c,
           yticks = [0,1],
           #xticksize = ticks,
           #yticksize = ticks,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\exp(-\Delta V)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Array of axes
axes = [ax2_1, ax2_2, ax2_3, ax2_4, ax2_5]

# Adjust whitespace between columns and rows
#colgap!(fig.layout, 120)
#rowgap!(fig.layout, 10)
