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
          #xticks = [-1,0,1],
          #yticks = [-0.35,0,0.5],
          xtickwidth = border,
          ytickwidth = border,
          xlabel = L"\textbf{parameter}",
          ylabel = L"\textbf{solution}",
          xlabelsize = labels,
          ylabelsize = labels,
          xticklabelsize = labels,
          yticklabelsize = labels,
         )
