include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data from csv
data = readin("../data/slide_03/May77.csv")
time = data[:,1]
parameter = data[:,2]
state = data[:,3]

printstyled("Generating the figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
# Create and customise the figure
fig, ax = mkfig(size = [1800,1200],
                pad = (30,80,10,30),
                box_position = [1,1],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((minimum(time),maximum(time)), (minimum(state),maximum(state))),
                title = L"\textbf{vegetation turbidity model}",
                toggle_title = true,
                title_color = RGBAf(0,0.236,0.236),
                lab = [L"\textbf{time}", L"\textbf{state}"],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                lab_pad = [-60.0,-50.0],
                ticks_lab_trunc = [0,1]
               )
# Customise the ticks for the plot
set_ticks(ax, time, state)
# Plot and export the figure 
lines!(ax, time, state, color = :teal, linewidth = 3.5)

save("../fig/slide_03_1.png", fig)
