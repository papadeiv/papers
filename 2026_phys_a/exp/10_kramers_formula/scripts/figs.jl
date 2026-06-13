"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 350))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Modified escape EWS 
ax1 = Axis(fig[1,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticksvisible = false,
           yticksvisible = false,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\exp(-\Delta V/D)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Prefactor 
ax2 = Axis(fig[1,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticksvisible = false,
           yticksvisible = false,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"(2\pi)^{-1}\sqrt{|V''(b)|V''(a)}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Kramer's escape formula 
ax3 = Axis(fig[1,3],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = false,
           ylabelvisible = false,
           xticksvisible = false,
           yticksvisible = false,
           xticklabelsvisible = false,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"R_{\text{a}\to\text{b}}",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Colorbar for the regularization coefficient at the bottom
Colorbar(fig[1,4],
         flipaxis = false,
         limits = (D_set[1], D_set[end]),
         spinewidth = border,
         ticksvisible = false,
         ticklabelsize = labels,
         ticks = [D_set[1], D_set[end]],
         tickwidth = border,
         tickalign = 1,
         labelsize = labels,
         colormap = :roma,
         vertical = true,
         width = 20.0,
         height = Relative(1.0)
        )

# Adjust whitespace between columns and rows
rowgap!(fig.layout, 20)
