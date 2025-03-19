include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

###########################
# End of the last ice age #
###########################

# Import the data from csv
data = readin("../data/slide_01/Petit01.csv")
time = data[:,1]
state = data[:,2]

# Create and customise the figure
fig, ax = mkfig(size = [1200,1200],
                title = L"\textbf{end of the last glaciation}",
                toggle_title = true,
                lab = [L"\textbf{kyears BCE}", L"\mathbf{^2}\textbf{H concentration}"],
                lab_pad = [-60.0,-60.0],
                ax_orientation = [true,false],
                x_ticks = [60,10],
                y_ticks = [-490,-425],
                ticks_lab_trunc = [0,0]
               )

# Plot and export the figure 
lines!(ax, time.*1e-3, state, color = :red, linewidth = 3.5)
save("../fig/1.png", fig)

#################################
# Desertification of the Sahara #
#################################

# Import the data from csv
data = readin("../data/slide_01/deMonecal00.csv")
time = data[:,1]
state = data[:,2]

# Create and customise the figure
fig, ax = mkfig(size = [1200,1200],
                title = L"\textbf{end of the African humid period}",
                toggle_title = true,
                lab = [L"\textbf{kyears BCE}", L"\textbf{sea surface temperature}"],
                lab_pad = [-60.0,-40.0],
                ax_orientation = [true,false],
                x_ticks = [14,2],
                y_ticks = [19,27],
                ticks_lab_trunc = [0,0]
               )

# Plot and export the figure 
lines!(ax, time.*1e-3, state, color = :green, linewidth = 3.5)
save("../fig/2.png", fig)

###############################
# Vegetation biomass collapse #
###############################

# Import the data from csv
data = readin("../data/slide_01/May77.csv")
time = data[:,1]
param = data[:,2]
state = data[:,3]

# Create and customise the figure
fig, ax = mkfig(size = [1200,1200],
                title = L"\textbf{vegetation biomass collapse}",
                toggle_title = true,
                lab = [L"\textbf{days}", L"\textbf{density}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [0,200],
                y_ticks = [0,8],
                ticks_lab_trunc = [0,0]
               )

# Plot and export the figure 
lines!(ax, LinRange(0.0,time[end],length(time)), state, color = :orange, linewidth = 3.5)
save("../fig/slide_01/3.png", fig)
