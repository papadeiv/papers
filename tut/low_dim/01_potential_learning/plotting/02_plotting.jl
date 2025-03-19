include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the data from csv
μ = readin("../data/parameter.csv")
constants = readin("../data/normalisation_constants.csv")
N = constants[:,1]
N_OUP = constants[:,2]

# Define the exact scalar potential function 
U(x, μ) = μ*x + x^2 - x^3 + (1/5)*x^4
# Define the exact stationary dprobability distribution
σ = 0.200::Float64
p(x, μ, N) = N*exp(-(2*U(x,μ))/(σ^2))

# Loop over the parameter values
printstyled("Generating figures of the least-squares solution\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        #########################
        # Approximate histogram #
        #########################

        # Import the data at the current parameter value
        distribution = readin("../data/distribution/$n.csv")

        # Define the limits for the equilibrium distribution histogram 
        x_inf = -(2.00::Float64)
        x_sup = 4.00::Float64
        y_inf = -(0.25::Float64)
        y_sup = 6.0::Float64 

        # Define the domain of the function
        domain = LinRange(x_inf, x_sup, 1000)

        # Create and customise the histogram figure 
        fig, ax = mkfig(size = [1200,900],
                        bg_out = "#eeeeeeff",
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\textbf{density}"],
                        lab_pad = [-60.0,-40.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [0,y_sup],
                        ticks_lab_trunc = [1,1]
                       )
        # Plot the stationary distribution
        lines!(ax, domain, [p(x,μ[n],N[n]) for x in domain], color = (:black,0.35), linewidth = 4.5)
        lines!(ax, domain, [p(x,μ[n],N_OUP[n]) for x in domain], color = (:blue,0.35), linewidth = 4.5)
        # Plot the histogram approximating the stationary distribution 
        scatter!(ax, distribution[:,1], distribution[:,2], color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the equilibrium distribution histogram
        save("../fig/distribution/$n.png", fig)

        #########################
        # Inverted distribution #
        #########################
         
        # Import the data at the current parameter value
        inverted_potential = readin("../data/potential/$n.csv")
        inverted_potential_OUP = readin("../data/potential/OUP$n.csv")
        coefficients = readin("../data/fit/$n.csv")

        # Define the least-squares polynomial fit 
        V = Polynomial(coefficients)

        # Define the limits for the scalar potential plot
        x_inf = -(1.75::Float64)
        x_sup = 4.75::Float64
        y_inf = -(0.5::Float64)
        y_sup = 3.5::Float64 

        # Define the domain of the function
        domain = LinRange(x_inf, x_sup, 1000)

        # Create and customise the scalar potential figure
        fig, ax = mkfig(size = [1200,900],
                        bg_out = "#eeeeeeff",
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        lab_pad = [-60.0,-60.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [y_inf,y_sup],
                        ticks_lab_trunc = [1,1]
                       )
        # Plot the scalar potential
        lines!(ax, domain, [U(x, μ[n]) for x in domain], color = (:black,0.25), linewidth = 4.5)
        # Plot the mean polynomial least-squares solution
        lines!(ax, domain, [V(x) for x in domain], color = (:red,1.00), linewidth = 4)
        # Plot the datapoints from the inverted distribution
        scatter!(ax, inverted_potential[:,1], inverted_potential[:,2], marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)
        scatter!(ax, inverted_potential_OUP[:,1], inverted_potential_OUP[:,2], marker = :xcross, color = (:blue,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the scalar potential function plot 
        save("../fig/fit/$n.png", fig)
end
