include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
equilibria = readin("../data/figure_07/equilibria.csv")
μ = equilibria[:,1]

# Import the solutions of the ensemble mean non-linear least-squares problem
c_nonlinear = readin("../data/figure_07/mean_coefficients.csv")
Nμ = length(c_nonlinear[:,1])

# Import the mean escape rates
escapes = readin("../data/figure_07/escape_rates.csv") 
escape_analytic = escapes[:,1]
escape_nonlinear = escapes[:,3]

# Import the escape ratios
escape_ratios = readin("../data/figure_07/escaped_percentage.csv")

# Import the variance ews 
variance_ews = readin("../data/figure_07/variance.csv")

# Import error estimates
estimates = readin("../data/figure_07/error_estimates.csv") 
error_nonlinear = estimates[:,2]

# Array to store the indices of the parameter values for which at least one trajectory in the ensemble survived
valid = Int64[]

# Loop over the parameter values
for n in 1:Nμ
        # Check if least one particle survived (didn't escape the basin)
        if escape_ratios[n] < 1.0::Float64
                # Get the current parameter index as valid
                push!(valid, n) 
        end
end

# Define index for the parameter value to start the plot
init = 100 

# Create and customise the error decay figure
fig, ax = mkfig(size = [1800,600],
                box_position = [1,1],
                limits = ((μ[init],μ[end]), (0,350)),
                lab = [L"\mathbf{\mu}", L"\mathbf{||V-V_{*}||_2}"],
                lab_size = [40,40],
                lab_pad = [-50.0,-60.0],
                x_ticks = [μ[init],μ[end]],
                y_ticks = [0,350],
                ticks_lab_size = [40,40],
                ticks_lab_xpos = [:right,:top],
                ticks_lab_trunc = [1,0]
)
lines!(ax, μ[valid], error_nonlinear[valid], color = (:teal,1.00), linewidth = 8)

# Create and customise the ensemble mean escape rate figure 
fig, ax = mkfig(fig = fig,
                box_position = [1,2],
                limits = ((μ[init],μ[end]), (-0.001,0.035)),
                lab = [L"\mathbf{\mu}", L"\mathbf{p_{esc}}"],
                lab_size = [40,40],
                lab_pad = [-50.0,-60.0],
                flip_y = true,
                x_ticks = [μ[init],μ[end]],
                y_ticks = [0, 0.035],
                ticks_lab_size = [40,40],
                ticks_lab_xpos = [:left,:top],
                ticks_lab_trunc = [1,3]
)
# Plot the escape rate for the non-linear reconstruction 
lines!(ax, μ[valid], escape_nonlinear[valid], color = (:teal,1.00), linewidth = 8)
# Plot the variance ews
lines!(ax, μ[valid], variance_ews[valid], color = (:darkgoldenrod2,0.75), linewidth = 8)

# Export the figure 
save("../fig/fig07.png", fig)


