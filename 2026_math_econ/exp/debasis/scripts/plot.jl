"""
    Plotting script
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the zero-sets of the two surfaces and their intersections (i.e. solutions of F(em_t+1, ew_t+1))
function plot_zero_set(Z)
        # Create empty layouts for the figures
        include("./scripts/figs.jl")

        # Extract the zero sets
        Z1, Z2 = Z

        # Plot the zero-sets
        contour!(ax, em_range, ew_range, Z1; levels = [0.0], color = CtpRed, linewidth = 4)
        contour!(ax, em_range, ew_range, Z2; levels = [0.0], color = CtpBlue, linewidth = 4)

        #=
        # Plot the intersections of the 0-sets
        scatter!(ax, -0.00411292, -0.0310015, markersize = 50, color = CtpMauve, strokewidth = 5.0, strokecolor = :black)
        text!(-0.00411292, -0.0310015, text = L"\mathbf{(-0.004, -0.031)}", color = CtpMauve, align = (:left, :bottom), fontsize = 30, offset = (30.0, 0.0))
        scatter!(ax, 3.24636, 2.88324, markersize = 50, color = CtpTeal, strokewidth = 5.0, strokecolor = :black)
        text!(3.24636, 2.88324, text = L"\mathbf{(3.246, \,2.883)}", color = CtpTeal, align = (:right, :baseline), fontsize = 30, offset = (-30.0, 0.0))
        =#

        savefig("$idx.png", fig)
end
