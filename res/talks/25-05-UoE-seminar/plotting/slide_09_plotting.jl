include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

solution = readin("../data/slide_09/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nμ = length(μ)

printstyled("Generating the figure\n"; bold=true, underline=true, color=:light_blue)
# Create and customise the potential axis 
fig, ax = mkfig(size = [1200,400],
                border = 0.0,
                limits = ((t[1], t[end]), (0.5, 1.25)),
                toggle_lab = [false, false],
                toggle_ticks = [false, false],
                toggle_ticks_lab = [false, false],
               )
# Plot the scalar potential 
lines!(ax, t, u, linewidth = 1, color = (:teal, 1.0))

# Export the figure
save("../fig/slide_09.png", fig)
