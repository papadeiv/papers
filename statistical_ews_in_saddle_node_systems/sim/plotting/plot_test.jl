include("../../inc/IO.jl")
include("../../inc/EscapeProblem.jl")
include("../../inc/SystemAnalysis.jl")
include("../../inc/PotentialLearning.jl")
include("../../inc/TimeseriesAnalysis.jl")

using CairoMakie, PrettyMakie, Makie.Colors
using LaTeXStrings

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

CtpMauve = colorant"rgb(202,158,230)"
CtpTeal = colorant"rgb(129, 200, 190)"
CtpBlue = colorant"rgb(140, 170, 238)"
CtpRed = colorant"rgb(231, 130, 132)"
CtpYellow = colorant"rgb(229,200,144)"

# Import the data from csv
coefficients = readin("../data/stationary_ensemble/coefficients.csv")
taylor = readin("../data/stationary_ensemble/taylor.csv")
solutions = readin("../data/stationary_ensemble/solutions.csv")
t = readin("../data/stationary_ensemble/time.csv")
Ne = length(solutions[:,1])
Nt = length(solutions[1,:])

# Define number of bins for the histograms of the random variables
Nb = convert(Int64, floor(0.02*Ne))

#=
##############
# Timeseries #
##############
        
# Create and customise the timeseries figure 
fig, ax = mkfig(size = [900,600],
                pad = (20,30,10,20), # Order is: left, right, bottom, top 
                bg_out = :white,
                box_position = [1,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((t[1], t[end]), (2.2, 2.8)),
                lab = [L"\mathbf{t}", L"\mathbf{x_t}"],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-50.0,-50.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [2.2,2.8],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [0,1]
               )
# Loop over the ensemble
@showprogress for n in 1:200
        # Extract the solution
        local u = solutions[n,:] 
        # Plot it
        lines!(ax, t, u, linewidth = 0.5, color = n, colormap = [(CtpMauve,0.25), (CtpTeal,0.25)], colorrange = (1,500))
end

# Export the timeseries figure
save("../fig/solutions.png", fig)
=#
 
################
# Coefficients #
################
        
# Assemble an histogram for c1
c1 = coefficients[:,1]
bins, pdf = fit_distribution(c1, n_bins=Nb+1)

# Create and customise the c1 distribution figure 
fig, ax = mkfig(size = [1800,600],
                pad = (20,75,10,20), # Order is: left, right, bottom, top 
                bg_out = :white,
                box_position = [1,1],
                border_color = colorant"rgb(76,79,105)",
                limits = ((bins[1], bins[end]), (0, nothing)),
                lab = [L"\mathbf{c_1}", L"\mathbf{p}"],
                toggle_lab = [true, false],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [bins[1],bins[end]],
                y_ticks = [0,maximum(pdf)],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [3,0]
               )
# Plot the histogram
barplot!(ax, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
# Plot the Taylor's coefficient
lines!(ax, [taylor[1], taylor[1]], [0, maximum(pdf)], linewidth = 3, color = :brown2)

# Assemble an histogram for c2
c2 = coefficients[:,2]
bins, pdf = fit_distribution(c2, n_bins=Nb+1)

# Create and customise the c2 distribution figure 
fig, ax = mkfig(fig = fig,
                box_position = [1,2],
                border_color = colorant"rgb(76,79,105)",
                limits = ((bins[1], bins[end]), (0, nothing)),
                lab = [L"\mathbf{c_2}", L"\mathbf{p(c_2)}"],
                toggle_lab = [true, false],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [bins[1],bins[end]],
                y_ticks = [0,maximum(pdf)],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [3,0]
               )
# Plot the histogram
barplot!(ax, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
# Plot the Taylor's coefficient
lines!(ax, [taylor[2], taylor[2]], [0, maximum(pdf)], linewidth = 3, color = :brown2)

# Assemble an histogram for c3
c3 = coefficients[:,3]
bins, pdf = fit_distribution(c3, n_bins=Nb+1)

# Create and customise the c3 distribution figure 
fig, ax = mkfig(fig = fig,
                box_position = [1,3],
                border_color = colorant"rgb(76,79,105)",
                limits = ((bins[1], bins[end]), (0, nothing)),
                lab = [L"\mathbf{c_3}", L"\mathbf{p(c_3)}"],
                toggle_lab = [true, false],
                lab_color = [colorant"rgb(76,79,105)", colorant"rgb(76,79,105)"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [bins[1],bins[end]],
                y_ticks = [0,maximum(pdf)],
                ticks_color = [colorant"rgb(76,79,105)",colorant"rgb(76,79,105)"],
                ticks_lab_trunc = [4,0]
               )
# Plot the histogram
barplot!(ax, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
# Plot the Taylor's coefficient
lines!(ax, [taylor[3], taylor[3]], [0, maximum(pdf)], linewidth = 3, color = :brown2)

# Export the coefficients' distribution figure
save("../fig/coefficients_distribution.png", fig)
