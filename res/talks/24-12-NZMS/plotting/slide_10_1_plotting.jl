include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
μ = readin("../data/slide_10_1/μ.csv")

# Define the limits for the scalar potential plot
x_inf = 2.0::Float64
x_sup = 4.0::Float64
y_inf = -(0.5::Float64)
y_sup = 6.0::Float64 

# Loop over the parameter values
printstyled("Generating figures of the evolution of the histrogram\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        # Import the data at the current parameter value
        local distribution = readin("../data/slide_10_1/distribution/$n.csv")

        # Create and customise the scalar potential figure
        fig, ax = mkfig(size = [1400,900],
                        bg_out = "#eeeeeeff",
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\textbf{density}"],
                        lab_pad = [-60.0,-40.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [0,y_sup],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the histogram approximating the stationary distribution 
        scatter!(ax, distribution[:,1], distribution[:,2], color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the scalar potential function plot 
        save("../fig/slide_10_1/$n.png", fig)
end
