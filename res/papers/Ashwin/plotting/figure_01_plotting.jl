include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

 # Import the data from csv 
solution = readin("../data/figure_01/solutions/1.csv")
t = solution[:,1]

# Create and customise the state axis 
fig, ax1 = mkfig(size = [1200,900],
                 box_position = [1,1],
                 limits = ((t[1], t[end]), (-2.1, 2.1)),
                )
 
# Create and customise the parameter axis 
fig, ax2 = mkfig(fig = fig,
                 box_position = [2,1],
                 limits = ((t[1], t[end]), (-1.10, 1.10)),
                )

# Define array of colors
colors = [:red, :blue, :green]

# Loop over the ramps
printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)
for n in 1:3
        # Import the data from csv 
        solution = readin("../data/figure_01/solutions/$n.csv")
        t = solution[:,1]
        μ = solution[:,2]
        u = solution[:,3]
        stabl_eq_1 = readin("../data/figure_01/stable_eq_1/$n.csv")
        stabl_eq_2 = readin("../data/figure_01/stable_eq_2/$n.csv")
        unstabl_eq = readin("../data/figure_01/unstable_eq/$n.csv")

        # Plot the bifurcation diagram of the frozen system
        if n == 1
                lines!(ax1, stabl_eq_1[:,1], stabl_eq_1[:,2], linewidth = 4, color = :black)
                lines!(ax1, stabl_eq_2[:,1], stabl_eq_2[:,2], linewidth = 4, color = :black)
                lines!(ax1, unstabl_eq[:,1], unstabl_eq[:,2], linewidth = 4, color = :black, linestyle = :dash)
        end
        # Plot the state's solution 
        lines!(ax1, t, u, linewidth = 6, color = (colors[n], 0.5))
        # Plot the parameter's shift
        lines!(ax2, t, μ, linewidth = 4, color = colors[n])

        # Export the figure 
        save("../fig/figure.png", fig)
end
