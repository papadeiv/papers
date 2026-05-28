"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

#=
# Specify the figure dimensions and style
fig = Figure(; size = (1200, 1200))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Variance EWS 
ax1 = Axis(fig[1,1], limits = (-1, 0.025, 0, 0.02),
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
ax2 = Axis(fig[2,1], limits = (-1, 0.025, 0.25, 1),
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
=#

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 800))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

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
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"\theta_1",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

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
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"\theta_2",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

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
             xtickwidth = border,
             ytickwidth = border,
             xlabel = L"\theta_3",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

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
             xlabel = L"\text{Var(x)}",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

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
             xlabel = L"\exp(-\Delta V)",
             ylabel = L"\text{distribution}",
             xlabelsize = labels,
             ylabelsize = labels,
             xticklabelsize = labels,
             yticklabelsize = labels,
            )

ax6 = Axis(fig[2,3], limits = (-2, 2, -1, 1),
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

# Array of axes
axes = [ax1, ax2, ax3, ax4, ax5, ax6]

# Colorbar for the regularization coefficient at the bottom
Colorbar(fig[3,1:3],
         limits = (α_set[1], α_set[end]),
         spinewidth = border,
         ticklabelsize = labels,
         ticks = α_set,
         tickwidth = border,
         label = L"\alpha",
         labelsize = labels,
         colormap = :phase,
         vertical = false,
         width = Relative(1.0),
         height = 20
        )

# Transparent grid overlaying the colorbar
colorbar = Axis(fig[3, 1:3],
                limits = (α_set[1], α_set[end], 0, 1),
                backgroundcolor = :transparent,
                xticksvisible = false,
                yticksvisible = false,
                xticklabelsvisible = false,
                yticklabelsvisible = false,
               )
hidedecorations!(colorbar)
hidespines!(colorbar)
