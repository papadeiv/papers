include("./mkfig.jl")

function mirror_axis(fig; 
                     box_position = [1,1],
                     color = :red,
                     toggle_lab = true,
                     ticks_lab_trunc = 1 
                    )
        nullfig, ax = mkfig(fig=fig,
                            box_position = box_position,
                            toggle_lab = [false, toggle_lab],
                            lab_color = [color, color],
                            flip_y = true,
                            ticks_lab_color = [color, color],
                            ticks_lab_trunc = [1, ticks_lab_trunc], 
                           )
        hidespines!(ax)
        hidexdecorations!(ax)
        return ax 
end
