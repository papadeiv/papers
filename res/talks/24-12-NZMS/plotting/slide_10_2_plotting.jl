include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
μ = readin("../data/slide_10_2/μ.csv")
c = readin("../data/slide_10_2/fit/mean_coefficients.csv")

# Define the exact scalar potential function 
U(x, μ) = μ*x + x^2 - x^3 + (1/5)*(x^4)

# Define a shifted analytical potential function for plotting purposes
Us(x, μ, x0) = U(x, μ) - U(x0, μ) 

# Define the limits for the scalar potential plot 
x_inf = -(2.0::Float64)
x_sup = 5.0::Float64
y_inf = -(4.0::Float64)
y_sup = 14.0::Float64

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Loop over the parameter values
printstyled("Generating figures of the least-squares solution\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        # Define the mean least-squares polynomial fit 
        V = Polynomial(c[n,:])

        # Get the stationary points of the interpolated polynomial
        points = get_stationary_points(V, domain=[3,4])

        #########################
        # Inverted distribution #
        #########################
         
        # Import the data at the current parameter value
        inverted_potential = readin("../data/slide_10_2/potential/$n.csv")

        # Create and customise the scalar potential figure
        fig, ax = mkfig(size = [1200,900],
                        bg_out = "#eeeeeeff",
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        lab_pad = [-60.0,-60.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [y_inf,y_sup],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the scalar potential
        lines!(ax, domain, [Us(x, μ[n], points[end]) for x in domain], color = (:black,1.0), linewidth = 4.5)
        # Plot the datapoints from the inverted distribution
        scatter!(ax, inverted_potential[:,1], inverted_potential[:,2], marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the scalar potential function plot 
        save("../fig/slide_10_2/potential/$n.png", fig)

        #########################
        # Polynomial regression #
        #########################
        
        # Create and customise the scalar potential figure
        fig, ax = mkfig(size = [1200,900],
                        bg_out = "#eeeeeeff",
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        lab_pad = [-60.0,-60.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [y_inf,y_sup],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the scalar potential
        lines!(ax, domain, [Us(x, μ[n], points[end]) for x in domain], color = (:black,0.25), linewidth = 4.5)
        # Plot the mean polynomial least-squares solution
        lines!(ax, domain, [V(x) for x in domain], color = (:red,1.00), linewidth = 4)
        # Plot the datapoints from the inverted distribution
        scatter!(ax, inverted_potential[:,1], inverted_potential[:,2], marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the scalar potential function plot 
        save("../fig/slide_10_2/fit/$n.png", fig)
end
