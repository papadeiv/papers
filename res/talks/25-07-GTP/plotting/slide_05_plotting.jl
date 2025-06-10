include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

printstyled("Generating the figures\n"; bold=true, underline=true, color=:darkgreen)

ews_var = colorant"rgb(136,57,239)"
ews_skw = colorant"rgb(230,69,83)"
CtpFlamingo = colorant"rgb(221, 120, 120)"
CtpYellow = colorant"rgb(229, 200, 144)"
CtpPeach = colorant"rgb(254, 100, 11)"
CtpTeal = colorant"rgb(23,146,153)"
CtpRed = colorant"rgb(210, 15, 57)"

##################
#   Timeseries   #
##################

# Import the data from csv
data = readin("../data/slide_05/solution.csv")
t = data[:,1] 
u = data[:,2]
Nt = length(t)

# Create and customise the timeseries figure
fig, ax = mkfig(size = [900,500],
                bg_out = colorant"rgb(239,241,245)",
                border_color = colorant"rgb(76,79,105)",
                limits = ((t[1], t[end]), (minimum(u), maximum(u))),
                lab = [L"\textbf{time}", L"\textbf{state}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the timeseries 
lines!(ax, t, u, color = CtpTeal, linewidth = 2.5)

# Export the timeseries figure 
save("../fig/slide_05.1.png", fig)

#################
#   Histogram   #
#################

# Import the data
histogram = readin("../data/slide_05/histogram.csv")
bins = histogram[:,1]
pdf = histogram[:,2]
Nb = length(bins)

# Create and customise the histogram figure
fig, ax = mkfig(size = [900,500],
                border_color = colorant"rgb(76,79,105)",
                limits = ((bins[1], bins[end]), (0, 10)),
                lab = [L"\textbf{state}", L"\mathbf{p}\textbf{(state)}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the histogram 
barplot!(ax, bins, pdf, color = pdf, colormap = [CtpYellow,CtpRed], strokecolor = :black, strokewidth = 1)

# Export the histogram figure 
save("../fig/slide_05.2.png", fig)

####################
#   Distribution   #
####################

# Define the noise level 
σ = 0.125::Float64
D = (σ^2)/2.0::Float64

# Import the data
coefficients = readin("../data/slide_05/coefficients.csv")
guess = coefficients[:,1]
c = coefficients[:,2]

# Define the domain of the distribution 
domain = LinRange(bins[1], bins[end], 1000)

# Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
Xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
c0(μ) = - μ[1]*Xs(μ) - μ[2]*(Xs(μ))^2 - μ[3]*(Xs(μ))^3

# Define an arbitrary cubic with the the above constraint on {c0}
V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)

# Define the stationary probability distribution
f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))
N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf), parameters=μ)
ρ(x, μ) = N(μ)*f(x, μ)

# Create and customise the distribution figure
fig, ax = mkfig(size = [900,500],
                border_color = colorant"rgb(76,79,105)",
                limits = ((bins[1], bins[end]), (0, 10)),
                lab = [L"\textbf{state}", L"\mathbf{p}\textbf{(state)}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the histogram 
barplot!(ax, bins, pdf, color = pdf, colormap = [colorant"rgb(239,241,245)", colorant"rgb(156,160,176)"], strokecolor = :black, strokewidth = 1)
# Plot the equilibrium solution on top of the previous figure 
lines!(ax, domain, [ρ(x, c) for x in domain], color = [ρ(x, c) for x in domain], colormap = [CtpYellow, CtpRed], linewidth = 8.5)

# Export the distribution figure 
save("../fig/slide_05.3.png", fig)

#################
#   Potential   #
#################

# Import the data
ys = readin("../data/slide_05/shift.csv")

# Define the limits for the scalar potential plot
x_inf = -(1.0::Float64)
x_sup = 3.5::Float64
y_inf = -(0.5::Float64)
y_sup = 2.5::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
U(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Define the fixed parameter value
μ0 = 0.625::Float64

# Define a shifted version of the reconstructed potential
Us(x, μ) = ys[1] + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)

# Create and customise the scalar potential figure
fig, ax = mkfig(size = [900,500],
                border_color = colorant"rgb(76,79,105)",
                limits = ((x_inf, x_sup), (y_inf, y_sup)),
                lab = [L"\textbf{state}", L"\mathbf{V}\textbf{(state)}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the true potential
lines!(ax, domain, [U(x,μ0) for x in domain], color = :teal, linewidth = 6.5)
# Plot the numerical reconstruction of the potential
lines!(ax, domain, [Us(x,c) for x in domain], color = CtpRed, linewidth = 6.5)

# Export the potential figure 
save("../fig/slide_05.4.png", fig)
