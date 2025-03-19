include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the data from csv
time = readin("../data/slide_06/time.csv")
u = readin("../data/slide_06/u.csv")

###################
# Full timeseries #
###################

# Create and customise the full timeseries figure
fig, ax = mkfig(size = [1200,600],
                bg_out = "#eeeeeeff",
                limits = ((0, time[end]), (-1, 3)),
                lab = [L"\textbf{time}", L"\textbf{state}"],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the state timeseries 
lines!(ax, LinRange(0.0, time[end], length(time)), u, linewidth = 3.5, color = (:black, 1.0))
 
# Export the full timeseries figure 
save("../fig/slide_06/1.png", fig)

############
# Subset 1 #
############

# Read the detrended timeseries (i.e. the residuals) 
ut = readin("../data/slide_06/ut.csv")

# Define the infimum and supremum indices for the subset
a = 1000::Int64
b = 3000::Int64

# Create and customise the subset timeseries subplot 
fig, ax1 = mkfig(size = [1200,400],
                 pad = 30,
                 bg_out = "#eeeeeeff",
                 box_position = [1,1:3],
                 limits = ((time[a], time[b]), (-0.15,0.15)),
                 toggle_lab = [false,false],
                 toggle_ticks = [false,false],
                 toggle_ticks_lab = [false,false],
                )
# Plot the subset timeseries 
lines!(ax1, LinRange(time[a], time[b], (b-a)+1), ut[a:b], linewidth = 1.0, color = (:black, 1.0))


# Create and customise the subset timeseries subplot 
nullfig, ax2 = mkfig(fig = fig,
                     box_position = [1,4],
                     limits = ((0,15), (-0.15,0.15)),
                     toggle_lab = [false,false],
                     toggle_ticks = [false,false],
                     toggle_ticks_lab = [false,false],
                    )
# Plot the histogram of the timeseries
hist!(ax2, ut[a:b], bins = 20, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)
# Export the subset 1 figure 
save("../fig/slide_06/2.png", fig)

############
# Subset 2 #
############

# Define the infimum and supremum indices for the subset
a = 14900::Int64
b = 16900::Int64

# Create and customise the subset timeseries subplot 
fig, ax1 = mkfig(size = [1200,400],
                 pad = 30,
                 bg_out = "#eeeeeeff",
                 box_position = [1,1:3],
                 limits = ((time[a], time[b]), (-0.15,0.15)),
                 toggle_lab = [false,false],
                 toggle_ticks = [false,false],
                 toggle_ticks_lab = [false,false],
                )
# Plot the subset timeseries 
lines!(ax1, LinRange(time[a], time[b], (b-a)+1), ut[a:b], linewidth = 1.0, color = (:black, 1.0))

# Create and customise the subset timeseries subplot 
nullfig, ax2 = mkfig(fig = fig,
                     box_position = [1,4],
                     limits = ((0,15), (-0.15,0.15)),
                     toggle_lab = [false,false],
                     toggle_ticks = [false,false],
                     toggle_ticks_lab = [false,false],
                    )
# Plot the histogram of the timeseries
hist!(ax2, ut[a:b], bins = 20, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)
# Export the subset 2 figure 
save("../fig/slide_06/3.png", fig)
