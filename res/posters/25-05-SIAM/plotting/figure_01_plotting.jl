include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

#################################
# Desertification of the Sahara #
#################################

# Import the data from csv
data = readin("../data/figure_01/deMonecal00.csv")
time = data[:,1]
state = data[:,2]

# Create and customise the figure
fig, ax = mkfig(size = [1800,600],
                box_position = [1,1],
                title = L"\textbf{end of the African humid period}",
                toggle_title = true,
                lab = [L"\textbf{kyears BCE}", L"\textbf{SST}"],
                lab_pad = [-60.0,-50.0],
                ax_orientation = [true,false],
                x_ticks = [14,2],
                y_ticks = [19,27],
                ticks_lab_trunc = [0,0]
               )

# Plot and export the figure 
lines!(ax, time.*1e-3, state, color = :teal, linewidth = 3.5)

################################
# 2008 global financial crisis #
################################

# Import the data from csv
Index = readin("../data/figure_01/Diks18_index.csv")
index = Index[150:end]

# Create and customise the figure
fig, ax = mkfig(fig=fig,
                box_position = [1,2],
                title = L"\textbf{2008-09 global financial crisis}",
                toggle_title = true,
                lab = [L"\textbf{days}", L"\textbf{S&P 500}"],
                lab_pad = [-60.0,-50.0],
                ax_orientation = [true,false],
                x_ticks = ([365,1], ["January '08", "June '09"]),
                y_ticks = [0.7,1.4],
               )

# Plot and export the figure 
lines!(ax, LinRange(1,length(index),length(index)), index.*1e-3, color = :darkgoldenrod2, linewidth = 3.5)
save("../fig/fig1.png", fig)
