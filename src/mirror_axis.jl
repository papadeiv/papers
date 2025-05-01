include("./mkfig.jl")

function mirror_axis(fig, x_lims; 
                     box_position = [1,1],
                     color = :red,
                     y_lab = L"\mathbf{y}",
                     toggle_lab = true,
                     lab_pad = -60.0,
                     ticks_lab_trunc = 2
                    )
        nullfig, mirror_ax = mkfig(fig=fig,
                                   box_position = box_position,
                                   limits = ((x_lims[1], x_lims[2]), nothing),
                                   toggle_lab = [false, toggle_lab],
                                   lab_color = [color, color],
                                   lab_pad =[0.0, lab_pad],
                                   flip_y = true,
                                   ticks_lab_color = [color, color],
                                   ticks_lab_trunc = [1, ticks_lab_trunc], 
                                  )
        mirror_ax.ylabel = y_lab
        hidespines!(mirror_ax)
        hidexdecorations!(mirror_ax)
        return mirror_ax 
end
