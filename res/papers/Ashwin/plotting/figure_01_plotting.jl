include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data from csv 
solution = readin("../data/figure_01/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]

stabl_eq_1 = readin("../data/figure_01/stable_eq_1.csv")
stabl_eq_2 = readin("../data/figure_01/stable_eq_2.csv")
unstabl_eq = readin("../data/figure_01/unstable_eq.csv")

# Create and customise the state axis 
fig, ax = mkfig(size = [1200,900],
                box_position = [1,1],
                limits = ((t[1], t[end]), (-2.1, 2.1)),
               )
# Plot the bifurcation diagram of the frozen system
lines!(ax, stabl_eq_1[:,1], stabl_eq_1[:,2], linewidth = 4, color = :black)
lines!(ax, stabl_eq_2[:,1], stabl_eq_2[:,2], linewidth = 4, color = :black)
lines!(ax, unstabl_eq[:,1], unstabl_eq[:,2], linewidth = 4, color = :black, linestyle = :dash)
# Plot the state's solution 
lines!(ax, t, u, linewidth = 6, color = (:red, 0.5))

# Create and customise the parameter axis 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                limits = ((t[1], t[end]), (-1.10, 1.10)),
               )
# Plot the parameter's shift
lines!(ax, t, μ, linewidth = 4, color = :black)

# Export the figure 
save("../fig/figure.png", fig)
