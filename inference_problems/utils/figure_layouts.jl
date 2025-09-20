"""
Utility functions to generate and customise complex figure layouts.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 29-08-2025
"""

#-----------------------#
#                       # 
#   figure_layouts.jl   #                     
#                       #
#-----------------------#


"""
$(TYPEDSIGNATURES)

Generate a 2-dimensional `Figure` and `Axis` struct given input parameters.

Now `Makie` is nice and all, but its default figure layouts are incredibly bad and definitely not publication quality.
As a result a lot of tweaking and customisation is necessary to generate each picture which gets kind of annoying to do each time.
With `makefig` we set a lot of the defaults to aesthetically pleasingly values while still providing ways to customising them.

## Keyword arguments
* `size=[1000,1000]`: half the number of pixels for each axis (so the final resolution is twice the amount in input) 
* `bg_out=:transparent`: color of the outer background of the figure 
* `pad=(60,60,30,30)`: padding of the plot from the `(left, right, bottom, top)` border
* `fig = Figure(; size = (size[1], size[2]), figure_padding = pad, backgroundcolor = bg_out)`: the figure object (useful when plotting subplots in the same figure)
* `box_position=[1,1]`: coordinate of the plot within the figure (useful when plotting subplots in the same figure)
* `bg_in=:white`: color of the background box of the plot
* `border=5.0`: thickness of the border (spine) 
* `border_color=:black`: color of the border (spine) 
* `limits=(nothing,nothing)`: limits of the x- and y-axis 
* `title=L"textbf{template title}"`: title at the top of the plot 
* `toggle_title=false`: toggles the visibility of the title 
* `title_size=50`: fontsize of the title 
* `title_color=:black`: color of the title 
* `title_gap=4.0`: distance between the title and the top border of the plot 
* `lab=[L"mathbf{x}",L"mathbf{y}"]`: labels of the x- and y-axis 
* `toggle_lab=[true,true]`: toggle the visibility of the labels 
* `lab_size=[50,50]`: fontsize of the labels 
* `lab_color=[:black,:black]`: color of the labels 
* `ax_scale=[identity,identity]`: scale of the x- and y-axis (can be `log10`)
* `lab_pad=[0.0,0.0]`: distance of the x- and y-axis label from the bottom and leftmost border of the x- and y-axis's ticks labels respectively 
* `ax_orientation=[false,false]`: orientation of the x- and y-axis; if `true` the axis will increase right to left for the x-axis an top to bottom for the y-axis 
* `flip_y=false`: position of the y-ax3s; if `true` the y-axes will be positioned on the right of the plot 
* `x_ticks=nothing`: ticks of the x-axes; defaults to `ax.xticks` 
* `y_ticks=nothing`: ticks of the y-axes; defaults to `ax.yticks` 
* `toggle_ticks=[true,true]`: toggles the visibility of the x- and y- axis ticks 
* `ticks_size=22`: length of the x- and y-axis ticks
* `ticks_color=[:black,:black]`: color of the x- and y-axis ticks 
* `ticks_thick=[5.0,5.0]`: thickness of the x- and y-axis ticks (should be the same unit as `border` for consistency) 
* `toggle_ticks_lab=[true,true]`: toggles the visibility if the x- and y-axis ticks labels 
* `ticks_lab_size=[50,50]`: fontsize of the x- and y-axis ticks labels (should be the same as `lab_size` for consistency) 
* `ticks_lab_xpos=[:center,:top]`: horizontal and vertical alignment of the ticks labels of the x-axis 
* `ticks_lab_color = [lab_color[1],lab_color[2]]`: color of the x- and y-axis ticks labels 
* `ticks_lab_trunc=[1,1]`: number of decimal digits displayed by the x- and y-axis ticks labels 

## Output
`fig, ax`

## Example
"""
function makefig(;
               size = [1000,1000],
               bg_out = :transparent,
               pad = (60,60,30,30), # Order is: left, right, bottom, top 
               fig = Figure(; size = (size[1], size[2]), figure_padding = pad, backgroundcolor = bg_out),
               box_position = [1,1],
               bg_in = :white,
               border = 5.0,
               border_color = :black,
               limits = (nothing,nothing),
               title = L"\textbf{template title}",
               toggle_title = false,
               title_size = 50,
               title_color = :black,
               title_gap = 4.0,
               lab = [L"\mathbf{x}",L"\mathbf{y}"],
               toggle_lab = [true,true],
               lab_size = [50,50],
               lab_color = [:black,:black],
               lab_pad = [0.0,0.0],
               ax_scale = [identity,identity],
               ax_orientation = [false,false],
               flip_y = false,
               x_ticks = nothing,
               y_ticks = nothing,
               toggle_ticks = [true,true],
               ticks_size = [22,22],
               ticks_color = [:black,:black],
               ticks_thick = [5.0,5.0],
               toggle_ticks_lab = [true,true],
               ticks_lab_size = [50,50],
               ticks_lab_xpos = [:center,:top],
               ticks_lab_color = [lab_color[1],lab_color[2]],
               ticks_lab_trunc = [1,1],
        )

        ax = Axis(fig[box_position[1],box_position[2]],
                # Background
                backgroundcolor = bg_in,
                spinewidth = border,
                leftspinecolor = border_color,
                topspinecolor = border_color,
                rightspinecolor = border_color,
                bottomspinecolor = border_color,
                xgridvisible = false,
                ygridvisible = false,
                limits = limits,
                # Title
                title = title,
                titlevisible = toggle_title,
                titlesize = title_size,
                titlecolor = title_color,
                titlealign = :center,
                titlegap = title_gap,
                # Axes labels
                xlabel = lab[1],
                ylabel = lab[2],
                xlabelvisible = toggle_lab[1],
                ylabelvisible = toggle_lab[2],
                xlabelsize = lab_size[1],
                ylabelsize = lab_size[2],
                xlabelcolor = lab_color[1],
                ylabelcolor = lab_color[2],
                xlabelpadding = lab_pad[1],
                ylabelpadding = lab_pad[2],
                # Axes scale, position and direction
                xscale = ax_scale[1], 
                yscale = ax_scale[2],
                xreversed = ax_orientation[1],
                yreversed = ax_orientation[2],
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xtickalign = 1,
                ytickalign = 1,
                xticksvisible = toggle_ticks[1],
                yticksvisible = toggle_ticks[2],
                xticksize = ticks_size[1],
                yticksize = ticks_size[2],
                xtickcolor = lab_color[1],
                ytickcolor = lab_color[2],
                xtickwidth = border,
                ytickwidth = border,
                # Ticks labels
                xticklabelsvisible = toggle_ticks_lab[1],
                yticklabelsvisible = toggle_ticks_lab[2],
                xticklabelsize = ticks_lab_size[1],
                yticklabelsize = ticks_lab_size[2],
                xticklabelalign = (ticks_lab_xpos[1],ticks_lab_xpos[2]),
                yticklabelalign = (:right,:center),
                xticklabelcolor = ticks_lab_color[1],
                yticklabelcolor = ticks_lab_color[2],
                xtickformat = "{:.$(ticks_lab_trunc[1])f}",
                ytickformat = "{:.$(ticks_lab_trunc[2])f}",
        )

        # Additional customisation
        if x_ticks != nothing
                ax.xticks = x_ticks
        end
        if y_ticks != nothing
                ax.yticks = y_ticks 
        end
        if flip_y == true
                ax.yaxisposition = :right
                ax.yticklabelalign = (:left,:center)
        end

        # Return the figure and axis handle
        return fig, ax
end

"""
$(TYPEDSIGNATURES)

Generate a `Figure` and `Axis` struct given input parameters.

Export a figure `fig` at the specified `path`.
"""
function savefig(path, fig)
        # Create the export directory if it doesn't exist
        fullpath = "../../res/fig/" * path 
        mkpath(dirname(fullpath))

        # Export the figure
        save(fullpath, fig)
end
