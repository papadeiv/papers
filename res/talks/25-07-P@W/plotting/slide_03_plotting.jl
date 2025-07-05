include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

CtpMauve = colorant"rgb(202,158,230)"
CtpRed = colorant"rgb(231,130,132)"

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

################
#   Solution   # 
################

# Import data from csv
solution = readin("../data/slide_03/solution.csv")
t = solution[:,1]
u = solution[:,3]

# Create and customise the solution figure
fig, ax = mkfig(size = [900,900],
                bg_out = colorant"rgb(239,241,245)",
                pad = (30,60,10,30), 
                border_color = colorant"rgb(76,79,105)",
                box_position = [1,1],
                limits = ((t[1], t[end]), (-0.8, 3.2)),
                lab = [L"\mathbf{t}",L"\mathbf{x_t}"],
                toggle_lab = [false,true],
                lab_pad = [-50.0,-80.0],
                lab_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                x_ticks = [t[1],t[end]],
                y_ticks = [-0.8,3.2],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,1],
               )
# Plot the timeseries of the solution 
lines!(ax, t, u, color = CtpMauve, linewidth = 1)

###########
#   EWS   # 
###########

# Import data from csv
EWS = readin("../data/slide_03/ews.csv")
t_ews = EWS[:,1]
u_ews = EWS[:,2]

# Create and customise the solution figure
fig, ax = mkfig(fig = fig,
                border_color = colorant"rgb(76,79,105)",
                box_position = [2,1],
                limits = ((t[1], t[end]), (0.017, 0.051)),
                lab = [L"\mathbf{t}",L"\textbf{Var} \mathbf{(x_t)}"],
                lab_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-70.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [0.017,0.051],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [0,1],
               )
# Plot the early warning signal timeseries
lines!(ax, t_ews, u_ews, color = CtpRed, linewidth = 3)

# Export the figure
save("../fig/slide_03.png", fig)
