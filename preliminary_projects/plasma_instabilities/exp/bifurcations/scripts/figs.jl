"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 1200))

# Set thickness and size of axes elements 
border = 2.0
labels = 20

# Timeseries
ax = Axis(fig[1,1], limits = (nothing, nothing, nothing, nothing), 
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
          xlabel = L"\textbf{time}",
          ylabel = L"\textbf{solution}",
          xlabelsize = labels,
          ylabelsize = labels,
          xticklabelsize = labels,
          yticklabelsize = labels,
         )

# Colorbar for the regularization coefficient at the bottom
Colorbar(fig[2,1],
         limits = (μ_set[1], μ_set[end]),
         spinewidth = border,
         ticklabelsize = labels,
         ticks = μ_set,
         tickwidth = border,
         label = L"\textbf{parameter}",
         labelsize = labels,
         colormap = :phase,
         vertical = false,
         width = Relative(1.0),
         height = 20
        )
