include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"

printstyled("Generating the figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)

################
#   Solution   # 
################

# Import the data from csv
solution = readin("../data/slide_04/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]

EWS = readin("../data/slide_04/ews.csv")
t_ews = EWS[1:Ns,1]

tipping_point = readin("../data/slide_04/tipping_point.csv")
Ns = convert(Int64, tipping_point[1]) - 1::Int64

# Create and customise the solution figure
fig, ax = mkfig(size = [1200,1200],
                bg_out = colorant"rgb(239,241,245)",
                pad = (30,60,10,30), 
                border_color = colorant"rgb(76,79,105)",
                box_position = [1,1],
                limits = ((t[1], t[end]), (-0.8, 3.2)),
                lab = [L"\mathbf{t}",L"\mathbf{x_t}"],
                toggle_lab = [false,true],
                lab_pad = [-50.0,10.0],
                lab_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                x_ticks = [t[1],t[end]],
                y_ticks = [-0.8,3.2],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,1],
               )
# Plot the location of the tipping point
lines!(ax, [t_ews[end], t_ews[end]], [-0.8, 3.2], linewidth = 5, linestyle = :dash, color = colorant"rgb(76,79,105)")
# Plot the timeseries of the solution 
lines!(ax, t, u, color = CtpMauve, linewidth = 3.5)

#############
#   Error   # 
#############

# Import the data from csv
error = readin("../data/slide_04/error.csv")

# Create and customise the error figure
fig, ax = mkfig(fig = fig,
                border_color = colorant"rgb(76,79,105)",
                box_position = [2,1],
                limits = ((t[1], t[end]), (0.1,10)),
                lab = [L"\mathbf{t}",L"\mathbf{||V-V_{*}||_{2}}"],
                toggle_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,0.0],
                ax_scale = [identity,log10],
                x_ticks = [t[1],t[end]],
                y_ticks = [0.1,1,10],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,1],
               )
# Plot the location of the tipping point
lines!(ax, [t_ews[end], t_ews[end]], [0.1, 10], linewidth = 5, linestyle = :dash, color = colorant"rgb(76,79,105)")
# Plot the reconstruction error timeseries
lines!(ax, t_ews, error[1:Ns], color = CtpTeal, linewidth = 3.5)

###########
#   EWS   # 
###########

# Import the data from csv
u_var = EWS[1:Ns,2]
u_esc = EWS[1:Ns,3]
true_esc = EWS[1:Ns,4]

# Create and customise the EWS figure
fig, ax = mkfig(fig = fig,
                border_color = colorant"rgb(76,79,105)",
                box_position = [3,1],
                limits = ((t[1], t[end]), (-1e-3, 0.016)),
                lab = [L"\mathbf{t}",L"\textbf{EWSs}"],
                lab_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,0.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [0, 0.016],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [0,2],
               )
# Plot the location of the tipping point
lines!(ax, [t_ews[end], t_ews[end]], [-1e-3, 0.016], linewidth = 5, linestyle = :dash, color = colorant"rgb(76,79,105)")
# Plot the early warning signals timeseries
lines!(ax, t_ews, true_esc, color = (:black,0.35), linewidth = 5)
lines!(ax, t_ews, u_var, color = CtpRed, linewidth = 4.5)
lines!(ax, t_ews, u_esc, color = CtpBlue, linewidth = 4.5)

# Export the figure
save("../fig/slide_04.png", fig)
