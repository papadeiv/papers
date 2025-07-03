include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

printstyled("Generating the figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)

#################################
# Desertification of the Sahara #
#################################

# Import the data from csv
data = readin("../data/figure_01/deMonecal00.csv")
time = data[:,1].*1e-3
state = data[:,2]

# Create and customise the figure
fig, ax = mkfig(size = [1800,600],
                pad = (30,80,10,30),
                box_position = [1,1],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((minimum(time),maximum(time)), (minimum(state),maximum(state))),
                title = L"\textbf{desertification of the Sahara}",
                toggle_title = true,
                title_color = RGBAf(0,0.236,0.236),
                lab = [L"\textbf{kyears BCE}", L"\textbf{SST}"],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                lab_pad = [-60.0,-50.0],
                ax_orientation = [true,false],
                ticks_lab_trunc = [0,0]
               )
# Customise the ticks for the plot
set_ticks(ax, time, state)
ax.xtickformat = "{:.0f}"
# Plot and export the figure 
lines!(ax, time, state, color = :teal, linewidth = 3.5)

################################
# 2008 global financial crisis #
################################

# Import the data from csv
Index = readin("../data/figure_01/Diks18_index.csv")
index = Index[150:end].*1e-3

# Create and customise the figure
fig, ax = mkfig(fig=fig,
                box_position = [1,2],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((1,365), (minimum(index),maximum(index))),
                title = L"\textbf{2008-09 financial crisis}",
                toggle_title = true,
                title_color = RGBAf(0,0.236,0.236),
                lab = [L"\textbf{days}", L"\textbf{S&P 500}"],
                lab_pad = [-60.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                ax_orientation = [true,false],
                x_ticks = ([365,1], ["Jan.'08", "Jan.'09"]),
                y_ticks = [minimum(index),maximum(index)],
               )

# Plot and export the figure 
lines!(ax, LinRange(1,length(index),length(index)), index, color = :teal, linewidth = 3.5)
save("../fig/slide_02.png", fig)
