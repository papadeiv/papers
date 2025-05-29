include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")

printstyled("Generating the figure for the escape early-warning signal\n"; bold=true, underline=true, color=:light_blue)

# Import the data from csv
residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])

solution = readin("../data/slides_10_to_13/solution.csv")
μ = solution[:,2]

escapes = readin("../data/slide_13/escapes.csv")
true_escape = escapes[:,1]
mean_escape = escapes[:,2]
lin_escape = escapes[:,3]

var_ews = readin("../data/slide_13/variance.csv")
μ_ews = var_ews[:,1]
var_ews = var_ews[:,2]

# Create and customise the escape rate in normal scale figure 
fig, ax = mkfig(size = [1200,800],
                box_position = [1,1],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((μ[1],μ[end]), (-0.01,0.3)),
                lab = [L"\mathbf{\mu}", L"\mathbf{p}_{\textbf{tip}}"],
                toggle_lab = [false,true],
                lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                lab_pad = [-50.0,-60.0],
                x_ticks = [μ[1],μ[Nw],μ[end]],
                y_ticks = [0,0.3],
                ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                toggle_ticks_lab = [false,true],
                ticks_lab_size = [40,40],
                ticks_lab_trunc = [1,4]
               )
# Plot the true escape rate 
lines!(ax, μ[Nw:end], true_escape, linewidth = 5, color = :black)
# Plot the variance ews
lines!(ax, μ_ews, var_ews, linewidth = 5, color = :orange)
# Plot the escape rate from the mean detrended timeseries 
lines!(ax, μ[Nw:end], mean_escape, linewidth = 5, color = colorant"rgb(207,19,92)")
# Plot the escape rate from the linear detrended timeseries 
lines!(ax, μ[Nw:end], lin_escape, linewidth = 5, color = colorant"rgb(38,106,176)")
 
# Create and customise the escape rate in log scale figure 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((μ[1],μ[end]), (1e-16,1.0)),
                lab = [L"\mathbf{\mu}", L"\mathbf{p}_{\textbf{tip}}"],
                lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                lab_pad = [-50.0,-60.0],
                ax_scale = [identity,log10],
                x_ticks = [μ[1],μ[Nw],μ[end]],
                y_ticks = [0,0.0003,0.3],
                ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                ticks_lab_size = [40,40],
                ticks_lab_trunc = [1,4]
               )
# Plot the true escape rate 
lines!(ax, μ[Nw:end], true_escape, linewidth = 5, color = :black)
# Plot the variance ews
lines!(ax, μ_ews, var_ews, linewidth = 5, color = :orange)
# Plot the escape rate from the mean detrended timeseries 
lines!(ax, μ[Nw:end], mean_escape, linewidth = 5, color = colorant"rgb(207,19,92)")
# Plot the escape rate from the linear detrended timeseries 
lines!(ax, μ[Nw:end], lin_escape, linewidth = 5, color = colorant"rgb(38,106,176)")
 
# Export the figure
save("../fig/slide_13.png", fig)
