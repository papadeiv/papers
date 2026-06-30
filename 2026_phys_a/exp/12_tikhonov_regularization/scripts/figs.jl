"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1200, 500))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Estimates in search space 
ax1 = Axis3(fig[1,1],
            aspect = (1,1,1),
            azimuth = 0.75*pi,
            viewmode = :fit,
            xspinewidth = border,
            yspinewidth = border,
            zspinewidth = border,
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
            xlabel = L"\theta_1^*",
            ylabel = L"\theta_2^*",
            zlabel = L"\theta_3^*",
            xlabelsize = labels,
            ylabelsize = labels,
            zlabelsize = labels,
            xticklabelsize = labels,
            yticklabelsize = labels,
            zticklabelsize = labels,
           )

# Regularization in search space 
ax2 = Axis3(fig[1,2],
            aspect = (1,1,1),
            azimuth = 0.575*pi,
            viewmode = :fit,
            xspinewidth = border,
            yspinewidth = border,
            zspinewidth = border,
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
            xlabel = L"\theta_1^*",
            ylabel = L"\theta_2^*",
            zlabel = L"\theta_3^*",
            xlabelsize = labels,
            ylabelsize = labels,
            zlabelsize = labels,
            xticklabelsize = labels,
            yticklabelsize = labels,
            zticklabelsize = labels,
           )

# Zoomed-in inset of the above search space 
ax3 = Axis3(fig[1,3],
            aspect = (1,1,1),
            azimuth = 0.85*pi,
            viewmode = :fit,
            xspinewidth = border,
            yspinewidth = border,
            zspinewidth = border,
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
            xlabel = L"\theta_1^*",
            ylabel = L"\theta_2^*",
            zlabel = L"\theta_3^*",
            xlabelsize = labels,
            ylabelsize = labels,
            zlabelsize = labels,
            xticklabelsize = labels,
            yticklabelsize = labels,
            zticklabelsize = labels,
           )

# Colorbar for the regularization coefficient at the bottom
Colorbar(fig[1,4],
         flipaxis = false,
         limits = (α_set[1], α_set[end]),
         spinewidth = border,
         ticksvisible = false,
         ticklabelsize = labels,
         ticks = [α_set[1], α_set[end]],
         tickwidth = border,
         tickalign = 1,
         labelsize = labels,
         colormap = :Accent_6,
         vertical = true,
         width = 20.0,
         height = Relative(0.625)
        )

# Adjust the figure separations between plots
colgap!(fig.layout, Relative(0.001))
