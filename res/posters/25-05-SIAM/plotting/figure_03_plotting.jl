include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

#############
# R-tipping #
#############

# Import the data from csv
time = readin("../data/figure_03/R-tipping/time.csv")
μ = readin("../data/figure_03/R-tipping/μ.csv")
u = readin("../data/figure_03/R-tipping/u.csv")

# Define the limits for the scalar potential plot
x_inf = -(25.0::Float64)
x_sup = 10.0::Float64
y_inf = -(1500.0::Float64)
y_sup = -y_inf 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = -x*(p^2 - 1.0) - p*x^2 - (1/3)*x^3

# Indices to export the figures
idx = [150, 480 ,492]

# Loop over the parameter values
printstyled("Generating figures for R-tipping\n"; bold=true, underline=true, color=:brown2)
@showprogress for n in 1:length(idx)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1000,1000],
                    border = 15.0,
                    pad = 30, # Equal amount of "whitespace" in all four directions 
                    limits = ((x_inf, x_sup), (y_inf, y_sup)),
                    toggle_lab = [false,false],
                    toggle_ticks = [false,false],
                    toggle_ticks_lab = [false,false],
                   )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[idx[n]]) for x in domain], color = :darkgoldenrod2, linewidth = 10)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[idx[n]], V(u[idx[n]],μ[idx[n]]), color = (:teal, 1.0), markersize = 75, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/fig3.1.$n.png", fig)
end

#############
# N-tipping #
#############

# Import the data from csv
time = readin("../data/figure_03/N-tipping/time.csv")
u = readin("../data/figure_03/N-tipping/u.csv")

# Define the limits for the scalar potential plot
x_inf = -(1.5::Float64)
x_sup = 3.5::Float64
y_inf = -(1.0::Float64)
y_sup = 3.5::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Define the fixed parameter value
μ0 = 1.0::Float64

# Indices to export the figures
idx = [20, 117, 340]
M = 10

# Loop over the parameter values
printstyled("Generating figures for N-tipping\n"; bold=true, underline=true, color=:darkgreen)
@showprogress for n in 1:length(idx)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1000,1000],
                    border = 15.0,
                    pad = 30, # Equal amount of "whitespace" in all four directions 
                    limits = ((x_inf, x_sup), (y_inf, y_sup)),
                    toggle_lab = [false,false],
                    toggle_ticks = [false,false],
                    toggle_ticks_lab = [false,false],
                   )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ0) for x in domain], color = :darkgoldenrod2, linewidth = 10)
        # Plot the previous states from the current timesteps
        for m in 1:M
                scatter!(ax, u[idx[n]-m], V(u[idx[n]-m],μ0), color = (:teal, 0.15), markersize = 70, strokewidth = 0)
        end
        # Plot the future states from the current timesteps
        for m in 1:M
                scatter!(ax, u[idx[n]+m], V(u[idx[n]+m],μ0), color = (:teal, 0.15), markersize = 70, strokewidth = 0)
        end
         # Plot the state (ball in the landscape)
        scatter!(ax, u[idx[n]], V(u[idx[n]],μ0), color = (:teal, 1.0), markersize = 70, strokewidth = 3)

        # Export the phase plane figure 
        save("../fig/fig3.2.$n.png", fig)
end
