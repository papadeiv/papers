include("../../../inc/IO.jl")
include("../../../inc/PlottingTools.jl")
include("../../../inc/EscapeProblem.jl")
include("../../../inc/SystemAnalysis.jl")
include("../../../inc/PotentialLearning.jl")
include("../../../inc/TimeseriesAnalysis.jl")

ews_var = colorant"rgb(136,57,239)"
ews_skw = colorant"rgb(230,69,83)"
CtpFlamingo = colorant"rgb(221, 120, 120)"
CtpYellow = colorant"rgb(229, 200, 144)"
CtpPeach = colorant"rgb(254, 100, 11)"
CtpTeal = colorant"rgb(23,146,153)"
CtpRed = colorant"rgb(210, 15, 57)"

# Import the data from csv
solution = readin("../data/sim_03/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(t)

residual = readin("../data/sim_03/residuals/1.csv")
Nw = length(residual)
w = trunc((Nw/Nt)*100, digits=1)

tipping_point = readin("../data/sim_03/tipping_point.csv")
Ns = convert(Int64, tipping_point[1])
tip_idx = convert(Int64, tipping_point[2])

histogram = readin("../data/sim_03/histograms/1.csv")
Nb = length(histogram[:,1])
b = trunc((Nb/Nw)*100, digits=1)

coefficients = readin("../data/sim_03/coefficients.csv")
Nc = length(histogram[1,:])

shifts = readin("../data/sim_03/shifts.csv")
shift_x = shifts[:,1]
shift_y = shifts[:,2]

EWS = readin("../data/sim_03/ews.csv")
μ_ews = EWS[1:Ns,1]
u_var = EWS[1:Ns,2]
u_esc = EWS[1:Ns,3]
true_esc = EWS[1:Ns,4]

# Define the noise level 
σ = 0.125::Float64
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
save("../fig/sim_03/ews.png", fig)

#=
####################
#   Distribution   #
####################

# Import the data from csv
histogram = readin("../data/sim_03/histograms/$Ns.csv")
bins = histogram[:,1]
pdf = histogram[:,2]

# Define plotting domain for the distribution
domain = collect(LinRange(bins[1], bins[end], 2000))
        
# Define the plotting limits for the potential
x_inf = -(1.5::Float64)
x_sup = 4.0::Float64
V_inf = -(1.0::Float64)
V_sup = 3.0::Float64

# Create and customise the distribution (and potential) figure 
fig, ax = mkfig(size = [600,600],
                border_color = colorant"rgb(76,79,105)",
                limits = ((x_inf, x_sup), (V_inf, V_sup)),
                lab = [L"\textbf{state}", L"\mathbf{V}\textbf{state}"],
                toggle_lab = [true,false],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [32.5,-50.0],
                flip_y = true,
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the stationary distribution of non-linear least-squares solution
lines!(ax, domain, [ρ(x, coefficients[end,:]) for x in domain], color = [ρ(x, coefficients[end,:]) for x in domain], colormap = [CtpYellow, CtpRed], linewidth = 5)

# Export the histogram fit figure
save("../fig/sim_03/distribution.png", fig)
 
#################
#   Potential   #
#################
        
# Define the shifted reconstructions
Vs(x) = shift_y[n] + coefficients[n,1]*(x - shift_x[n]) + coefficients[n,2]*(x - shift_x[n])^2 + coefficients[n,3]*(x - shift_x[n])^3

# Define plotting domain for the histogram
domain = collect(LinRange(x_inf, x_sup, 2000))
# Plot the true potential
lines!(ax, domain, [U(x, μ[n_end]) for x in domain], color = :teal, linewidth = 6.5)
# Plot the true potential
lines!(ax, domain, [Vs(x) for x in domain], color = CtpRed, linewidth = 6.5)

# Export the scalar potential figure
save("../fig/sim_03/potential.png", fig)
=#
