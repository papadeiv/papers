include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

CtpMauve = colorant"rgb(202,158,230)"
CtpYellow = colorant"rgb(229,200,144)"

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

# Define the limits for the scalar potential plot
x_inf = -(1.0::Float64)
x_sup = 3.5::Float64
y_inf = -(1.0::Float64)
y_sup = 3.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

################
#   Figure 1   # 
################

# Create and customise the potential and distribution figure
fig, ax = mkfig(size = [600,600],
                bg_out = colorant"rgb(239,241,245)",
                border_color = colorant"rgb(76,79,105)",
                pad = 30, # Equal amount of "whitespace" in all four directions 
                limits = ((x_inf, x_sup), (y_inf, y_sup)),
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )

# Import data from csv
solution = readin("../data/slide_02/solution.csv")
t = solution[:,1]
u = solution[:,2]

# Define the level of additive noise
σ = 0.50::Float64
D = (σ^2)/2

# Define the bifurcation parameter value
μ = 0.6::Float64

# Define the scalar potential 
V(x) = μ*x + x^2 - x^3 + (1/5)*(x^4)
 
# Define the exact stationary probability distribution
f(x, q) = exp(-(1.0::Float64/D)*(V(x)))
N = get_normalisation_constant(f, (-Inf, Inf))
ρ(x, q) = N*f(x, q)

# Plot the scalar potential
lines!(ax, domain, [V(x) for x in domain], color = :black, linewidth = 6)
# Plot the ball in the potential
scatter!(ax, u[end], V(u[end]), color = (CtpMauve, 1.0), markersize = 45, strokewidth = 3)

# Export figure 1 
save("../fig/slide_02.1.png", fig)

################
#   Figure 2   # 
################

# Plot the density function 
lines!(ax, domain, [ρ(x, 0) for x in domain], color = [ρ(x, 0) for x in domain], colormap = [CtpYellow, CtpMauve], linewidth = 6)

# Export figure 2 
save("../fig/slide_02.2.png", fig)

################
#   Figure 3   # 
################

# Import the data from csv
histogram = readin("../data/slide_02/histogram.csv")
bins = histogram[:,1]
pdf = histogram[:,2]

# Plot the histogram at the endtime 
barplot!(ax, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
# Replot the objects in Figures 1 and 2 so that they appear on top of the histogram
lines!(ax, domain, [V(x) for x in domain], color = :black, linewidth = 6)
scatter!(ax, u[end], V(u[end]), color = (CtpMauve, 1.0), markersize = 45, strokewidth = 3)
lines!(ax, domain, [ρ(x, 0) for x in domain], color = [ρ(x, 0) for x in domain], colormap = [CtpYellow, CtpMauve], linewidth = 6)

# Export figure 3 
save("../fig/slide_02.3.png", fig)

################
#   Figure 4   # 
################

# Define new plotting domain for the timeseries
new_domain = LinRange(x_inf, x_sup, Nt)

# Create and customise the timeseries figure
fig, ax = mkfig(size = [600,300],
                bg_out = colorant"rgb(239,241,245)",
                border_color = colorant"rgb(76,79,105)",
                pad = 30, # Equal amount of "whitespace" in all four directions 
                limits = ((t[1], t[end]), (x_inf, x_sup)),
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )

# Plot the timeseries
lines!(ax, t, u, color = [ρ(x, 0) for x in u], colormap = [CtpYellow, CtpMauve], linewidth = 4)

# Export figure 4
save("../fig/slide_02.4.png", fig)
