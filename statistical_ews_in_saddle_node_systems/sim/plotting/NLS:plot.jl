include("../../inc/IO.jl")
include("../../inc/PlottingTools.jl")
include("../../inc/EscapeProblem.jl")
include("../../inc/SystemAnalysis.jl")
include("../../inc/PotentialLearning.jl")
include("../../inc/TimeseriesAnalysis.jl")

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

ews_var = colorant"rgb(136,57,239)"
ews_skw = colorant"rgb(230,69,83)"
CtpFlamingo = colorant"rgb(221, 120, 120)"
CtpYellow = colorant"rgb(229, 200, 144)"
CtpPeach = colorant"rgb(254, 100, 11)"
CtpTeal = colorant"rgb(23,146,153)"
CtpRed = colorant"rgb(210, 15, 57)"

# Import the data from csv
solution = readin("../data/NLS/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(t)

residual = readin("../data/NLS/residuals/1.csv")
Nw = length(residual)
w = trunc((Nw/Nt)*100, digits=1)

tipping_point = readin("../data/NLS/tipping_point.csv")
tip_idx = convert(Int64, tipping_point[1])

histogram = readin("../data/NLS/histograms/1.csv")
Nb = length(histogram[:,1])
b = trunc((Nb/Nw)*100, digits=1)

coefficients = readin("../data/NLS/coefficients.csv")
Nc = length(histogram[1,:])

shifts = readin("../data/NLS/shifts.csv")
shift_x = shifts[:,1]
shift_y = shifts[:,2]

EWS = readin("../data/NLS/ews.csv")
μ_ews = EWS[:,1]
u_var = EWS[:,2]
u_esc = EWS[:,3]
true_esc = EWS[:,4]

# Define the noise level 
σ = 0.250::Float64
D = (σ^2)/2.0::Float64

# Define the true potential
U(x, μ) = x*μ + x^2 - x^3 + (1/5)*(x^4)

# Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
Xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
c0(μ) = - μ[1]*Xs(μ) - μ[2]*(Xs(μ))^2 - μ[3]*(Xs(μ))^3

# Define an arbitrary cubic potential with the the above constraint on {c0}
V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)

# Define the stationary probability distribution
f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))
N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf), parameters=μ)
ρ(x, μ) = N(μ)*f(x, μ)

##############
# Timeseries #
##############
        
# Define plotting limits for the timeseries
x_inf = -(1.0::Float64) 
x_sup = 3.0::Float64 

# Create and customise the timeseries figure 
fig, ax = mkfig(size = [1500,1000],
                pad = (50,50,10,10), # Order is: left, right, bottom, top 
                bg_out = :white,
                box_position = [1,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (x_inf, x_sup)),
                title = L"\mathbf{N_t = %$Nt}\,\textbf{, }\;\mathbf{N_w = %$Nw}\;\mathbf{(%$w %)}\,\textbf{, }\;\mathbf{N_b = %$Nb}\;\mathbf{(%$b %)}",
                toggle_title = true,
                title_color = colorant"rgb(76,79,105)",
                title_gap = 14.0,
                lab = [L"\mathbf{\mu}", L"\mathbf{x_t}"],
                toggle_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [x_inf,x_sup],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,1]
               )
# Plot the stationary timeseries
lines!(ax, μ, u, linewidth = 2.5, color = :teal)
# Plot the sliding window
poly!(ax, Point2f[(μ[1], x_inf), (μ[Nw], x_inf), (μ[Nw], x_sup), (μ[1], x_sup)], color = (colorant"rgb(147,237,213)", 0.35), strokecolor = :grey, strokewidth = 0.05)
# Plot the location of the tipping point 
lines!(ax, [μ[tip_idx], μ[tip_idx]], [x_inf, x_sup], linewidth = 5, linestyle = :dash, color = (:black, 0.4))

############
#   EWSs   #
############

# Define plotting limits for the EWSs 
y_inf = -5e-4 
y_sup = 0.012::Float64 

# Create and customise the EWSs figure 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (y_inf, y_sup)),
                lab = [L"\mathbf{\mu}", L"\textbf{signal}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-92.5],
                x_ticks = [μ[1],μ[tip_idx],μ[end]],
                y_ticks = [0,y_sup],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [2,3]
               )
# Plot the sliding window
poly!(ax, Point2f[(μ[1], y_inf), (μ[Nw], y_inf), (μ[Nw], y_sup), (μ[1], y_sup)], color = (colorant"rgb(147,237,213)", 0.35), strokecolor = :grey, strokewidth = 0.05)
# Plot the location of the tipping point 
lines!(ax, [μ[tip_idx], μ[tip_idx]], [y_inf, y_sup], linewidth = 5, linestyle = :dash, color = (:black, 0.4))
# Plot the EWSs 
lines!(ax, μ_ews, u_var, linewidth = 3.5, color = ews_var)
lines!(ax, μ_ews, true_esc, linewidth = 6, color = (:teal,0.5))
lines!(ax, μ_ews, u_esc, linewidth = 3.5, color = ews_skw)

# Export the EWS figure
save("../fig/NLS/ews.png", fig)

#############
#   Error   #
#############

# Import the data from csv
error = readin("../data/NLS/error.csv")
taylor = readin("../data/NLS/taylor.csv")

# Create and customise the error figure 
fig, ax = mkfig(size = [2400,1200],
                bg_out = :white,
                box_position = [1:4,1:3],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (0, 7)),
                lab = [L"\mathbf{\mu}", L"\mathbf{||V-V_{*}||_2}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [0, 7]
               )
# Plot the numerical error 
lines!(ax, μ[Nw:tip_idx], error, color = :black, linewidth = 5)

# Create and customise the 1-st Taylor coefficient subfigure
fig, ax = mkfig(fig = fig,
                box_position = [5:6,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (minimum(coefficients[:,1]), maximum(coefficients[:,1]))),
                lab = [L"\mathbf{\mu}", L"\mathbf{c_{1}}"],
                toggle_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                toggle_ticks = [false,true],
                toggle_ticks_lab = [false,true],
                y_ticks = [minimum(coefficients[:,1]), maximum(coefficients[:,1])]
               )
# Plot the Taylor coefficients timeseries 
lines!(ax, μ[Nw:tip_idx], taylor[:,2], color = :black, linewidth = 5)
# Plot the optimised coefficient solution 
lines!(ax, μ[Nw:tip_idx], coefficients[:,1], color = :red, linewidth = 5)

# Create and customise the 2-nd Taylor coefficient subfigure
fig, ax = mkfig(fig = fig,
                box_position = [5:6,2],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (minimum(coefficients[:,2]), maximum(coefficients[:,2]))),
                lab = [L"\mathbf{\mu}", L"\mathbf{c_{2}}"],
                toggle_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                toggle_ticks = [false,true],
                toggle_ticks_lab = [false,true],
                y_ticks = [minimum(coefficients[:,2]), maximum(coefficients[:,2])]
               )
# Plot the Taylor coefficients timeseries 
lines!(ax, μ[Nw:tip_idx], taylor[:,3], color = :black, linewidth = 5)
# Plot the optimised coefficient solution 
lines!(ax, μ[Nw:tip_idx], coefficients[:,2], color = :red, linewidth = 5)

# Create and customise the 3-rd Taylor coefficient subfigure
fig, ax = mkfig(fig = fig,
                box_position = [5:6,3],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (minimum(coefficients[:,3]), maximum(coefficients[:,3]))),
                lab = [L"\mathbf{\mu}", L"\mathbf{c_{3}}"],
                toggle_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                toggle_ticks = [false,true],
                toggle_ticks_lab = [false,true],
                y_ticks = [minimum(coefficients[:,3]), maximum(coefficients[:,3])]
               )
# Plot the Taylor coefficients timeseries 
lines!(ax, μ[Nw:tip_idx], taylor[:,4], color = :black, linewidth = 5)
# Plot the optimised coefficient solution 
lines!(ax, μ[Nw:tip_idx], coefficients[:,3], color = :red, linewidth = 5)

# Export the error figure
save("../fig/NLS/error.png", fig)

################
#   Analysis   #
################

# Import the data from csv
parameters = readin("../data/NLS/parameters.csv")

# Create and customise the second derivative figure 
fig, ax = mkfig(size = [1800,1600],
                bg_out = :white,
                box_position = [1,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (minimum(parameters[:,6]), maximum(parameters[:,5]))),
                lab = [L"\mathbf{\mu}", L"\mathbf{V^{''}}"],
                toggle_lab = [false,true],
                toggle_ticks_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [minimum(parameters[:,6]), maximum(parameters[:,5])]
               )
# Plot the second derivatives of the ground truth 
lines!(ax, μ[Nw:tip_idx], parameters[:,1], color = :black, linewidth = 5)
lines!(ax, μ[Nw:tip_idx], parameters[:,2], color = :black, linewidth = 5)
# Plot the second derivatives of the approximation
lines!(ax, μ[Nw:tip_idx], parameters[:,5], color = :red, linewidth = 5)
lines!(ax, μ[Nw:tip_idx], parameters[:,6], color = :blue, linewidth = 5)

# Create and customise the potential barrier figure 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (-0.5, maximum(parameters[:,7]))),
                lab = [L"\mathbf{\mu}", L"\mathbf{\Delta V}"],
                toggle_lab = [false,true],
                toggle_ticks_lab = [false,true],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [-0.5, maximum(parameters[:,7])]
               )
# Plot the barrier of the ground truth 
lines!(ax, μ[Nw:tip_idx], parameters[:,3], color = :black, linewidth = 5)
# Plot the barrier of the approximation
lines!(ax, μ[Nw:tip_idx], parameters[:,7], color = :red, linewidth = 5)

# Create and customise the escape ews figure 
fig, ax = mkfig(fig = fig,
                box_position = [3,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((μ[1], μ[end]), (-1e-3, 0.08)),
                lab = [L"\mathbf{\mu}", L"\textbf{exp}\mathbf{(-2\Delta V/\sigma^2)}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [0, 0.08]
               )
# Plot the escape rate of the ground truth 
lines!(ax, μ[Nw:tip_idx], parameters[:,4], color = :black, linewidth = 5)
# Plot the escape rate of the approximation
lines!(ax, μ[Nw:tip_idx], parameters[:,8], color = :red, linewidth = 5)

# Export the analysis figure
save("../fig/NLS/analysis.png", fig)
