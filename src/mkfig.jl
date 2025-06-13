# Create a customised figure for scientific plots 

CairoMakie.activate!(; px_per_unit = 2)

function mkfig(;
               size = [1000,1000],
               alpha = 1.0,
               bg_out = :transparent,
               pad = (30,60,10,30), # Order is: left, right, bottom, top 
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
               lab_color = [(:black,alpha),(:black,alpha)],
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
