include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

printstyled("Generating the figure\n"; bold=true, underline=true, color=:light_blue)

 # Import the data from csv 
solution = readin("../data/figure_03/solution.csv")
t = solution[:,1]
λ = solution[:,2]

stable_eq_1 = readin("../data/figure_03/stable_eq_1.csv")
stable_eq_2 = readin("../data/figure_03/stable_eq_2.csv")
unstable_eq = readin("../data/figure_03/unstable_eq.csv")

# Create and customise the state axis 
fig, ax = mkfig(size = [1200,900],
                pad = (60,40,10,30), # Order is: left, right, bottom, top 
                bg_out = :white,
                box_position = [1,1],
                limits = ((t[1], t[end]), (-2.1, 2.1)),
                lab = [L"\mathbf{t}",L"\mathbf{x\,(\lambda(t))}"],
                lab_size = [30,30],
                lab_pad = [-30.0,-30.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [-2,2],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Plot the bifurcation diagram of the frozen system
lines!(ax, stable_eq_1[:,1], stable_eq_1[:,2], linewidth = 4, color = :black)
lines!(ax, stable_eq_2[1:403,1], stable_eq_2[1:403,2], linewidth = 4, color = :black)
lines!(ax, stable_eq_2[404:801,1], stable_eq_2[404:801,2], linewidth = 4, color = :black)
lines!(ax, stable_eq_2[802:end,1], stable_eq_2[802:end,2], linewidth = 4, color = :black)
lines!(ax, unstable_eq[1:403,1], unstable_eq[1:403,2], linewidth = 4, color = (:black,0.35))
lines!(ax, unstable_eq[404:801,1], unstable_eq[404:801,2], linewidth = 4, color = (:black,0.35))
lines!(ax, unstable_eq[802:end,1], unstable_eq[802:end,2], linewidth = 4, color = (:black,0.35))
# Plot the state's solution 
lines!(ax, t, u, linewidth = 3, color = (:brown2, 0.75))

# Create and customise the parameter axis 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                limits = ((t[1], t[end]), (-2.10, 2.10)),
                lab = [L"\mathbf{t}",L"\mathbf{\lambda(t)\in\mathcal{P}(-2,+2)}"],
                lab_size = [30,30],
                lab_pad = [-30.0,-30.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [-2.0,2.0],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Plot the parameter's shift
lines!(ax, t, μ, linewidth = 4, color = :black)

# Export the figure 
save("../fig/figure_03.png", fig)
