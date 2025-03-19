include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data from csv
index = readin("../data/abstract/Diks18_index.csv")

# Create and customise the figure
fig, ax = mkfig(size = [1200,1200],
                title = L"\textbf{2008-09 global financial crisis}",
                toggle_title = true,
                lab = [L"\textbf{days}", L"\textbf{S&P 500 (kUSD)}"],
                lab_pad = [-60.0,-40.0],
                ax_orientation = [true,false],
                x_ticks = ([495,20], ["February '08", "January '10"]),
                y_ticks = [0.7,1.4],
               )

# Plot and export the figure 
lines!(ax, LinRange(1,length(index),length(index)), index.*1e-3, color = :green, linewidth = 3.5)
save("../fig/abstract/fig2.png", fig)
