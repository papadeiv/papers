include("../../../../../inc/PlottingTools.jl")

# Define the parametric family of solutions of the ODE
y(x, C) = (4.0/5.0)*x^2 - (3.0/4.0)*x + C/(x^3)

# Define the domain of the solution
x_inf = 0.05::Float64 
x_sup = 1.00::Float64
domain = LinRange(x_inf,x_sup,1000)

# Define a range of values for the integration constant
C_inf = -(1.0::Float64)
C_sup = 1.0::Float64
C = LinRange(C_inf,C_sup,100)

# Create and customise the figure
y_inf = -(2.0::Float64)
y_sup = 2.0::Float64
fig, ax = mkfig(size = [1800,1200],
                limits = ((x_inf,x_sup+0.009), (y_inf,y_sup)),
                lab = [L"\mathbf{x}", L"\mathbf{y}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [x_inf,x_sup],
                y_ticks = [y_inf,y_sup],
                ticks_lab_trunc = [1,1]
               )
# Loop over the constant values
for c in C
        # Plot the solution at the current constant value
        lines!(ax, domain, [y(x,c) for x in domain], color = (:black,0.25), linewidth = 4)
end
# Plot the only bounded solution 
lines!(ax, domain, [y(x,(0.0)) for x in domain], color = (:darkgreen,1.0), linewidth = 4)
# Plot the solution of the IVP
lines!(ax, domain, [y(x,(-1.0/20.0)) for x in domain], color = (:brown2,1.0), linewidth = 6)
# Plot the BC
scatter!(ax, 1.0, 0.0, color = (:brown2,1.0), strokecolor = :black, strokewidth = 1, markersize = 35)

# Export the figure
save("../fig/figure1.2.0.0.png", fig)
