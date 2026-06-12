"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 800))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Marozov's discrepancy principle
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
          xlabel = L"\alpha",
          ylabel = L"\text{misfit}",
          xlabelsize = labels,
          ylabelsize = labels,
          xticklabelsize = labels,
          yticklabelsize = labels,
         )
