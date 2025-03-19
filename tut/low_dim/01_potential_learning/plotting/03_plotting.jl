include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/IO.jl")


# Define the analytical potential function
U(x, μ) = -(1.0/3.0)*x^3 - μ*x

# Define a shifted analytical potential function for plotting purposes
Us(x, μ, x0) = U(x, μ) - U(x0, μ) 

# Import the parameter range
μ = readin("../data/parameter.csv")

# Import the ensemble mean and variance of the least-squares coefficients
coeff_mean = readin("../data/coefficients/mean.csv")
coeff_var = readin("../data/coefficients/variance.csv")

# Get the number of parameter's values
Nμ = length(μ)-1

#################################################
#       Potential function from data 
#################################################

# Loop over the parameter values
printstyled("Creating the plots\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Read in the histograms and inverted potential function of one reppresentative trajectory
        distribution = readin("../data/distribution/fixed_monostable_saddle_node_$n.csv")
        bins = distribution[:,1]
        pdf = distribution[:,2]
        potential = readin("../data/potential/fixed_monostable_saddle_node_$n.csv")
        xs = potential[:,1]
        Vs = potential[:,2]

        # Read in the settings of the simulation
        settings = readin("../data/settings/$n.csv")
        local σ = settings[1]
        local m = settings[2]
        local a = settings[3]
        local b = settings[4]

        # Define the left edges of the bins of the histogram
        I = b - a
        Nbins = length(bins)
        x = LinRange(a - 0.05*I, b + 0.05*I, Nbins+1)

        # Define the ploynomial least-square regression of the potential function
        V = Polynomial(coeff_mean[n,:])

        # Find the stationary points on the polynomial
        points = get_stationary_points(V)

        # Define the range of the potential function for plotting purposes
        A = min(a, points[1])
        B = max(b, points[end])

        # Create the figure and the subplot axes
        local fig, ax = rowfig(2, 1, 
                         size = (2000,1000),
                         x_ticks_lab_trunc = [2,2],
                         y_ticks_lab_trunc = 2
                        )

        # Customise the ticks in the subplot
        set_ticks(ax[1], bins, pdf)
        set_ticks(ax[2], xs, Vs)

        # Plot the equilibrium gaussian distribution
        lines!(ax[1], LinRange(a - 0.05*I, b + 0.05*I, 1000), [gaussian(x, m, sqrt(σ)) for x in LinRange(a - 0.05*I, b + 0.05*I, 1000)], color = (:black, 0.5), linewidth = 3)
        # Plot the empirical distribution
        scatter!(ax[1], bins, pdf, color = (:gray69,0.01), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Plot the actual (shifted) potential function 
        lines!(ax[2], LinRange(A - 0.05*I, B + 0.05*I, 1000), [Us(x, μ[n], m) for x in LinRange(A - 0.05*I, B + 0.05*I, 1000)], color = (:black, 0.5), linewidth = 3)
        # Plot the reconstructed potential function through least-square regression
        lines!(ax[2], LinRange(A - 0.05*I, B + 0.05*I, 1000), [V(x) for x in LinRange(A - 0.05*I, B + 0.05*I, 1000)], color = (:black, 1.0), linewidth = 3, linestyle = :dash)
        # Plot the inverted potential function
        scatter!(ax[2], xs, Vs, marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)
        #Plot the bins edges on both subplots
        scatter!(ax[1], x, zeros(Nbins+1), markersize = 5, color = :red)
        scatter!(ax[2], x, zeros(Nbins+1), markersize = 5, color = :red)

        # Export the figure
        save("../fig/potential/$n.png", fig)
end

#################################################
#       Coefficients mean and variance
#################################################

# Get the number of coefficients
n_coeff = size(coeff_mean,2)

# Create the figure and the subplot axes for the mean
fig, mean_axs = colfig(n_coeff, 1, size = [2000,2000], x_ticks_lab_trunc = 1)
# Create the subplot axes for the variance
fig, var_axs = colfig(n_coeff, 2, fig = fig, flip_y = true, x_ticks_lab_trunc = 1)

# Define scaling for the y-axes
mean_scale_value = [1e0,1e0,1e0,1e0] 
mean_scale_label= [L"\mathbf{\times 10^{5}}", "", "", ""] 
var_scale_value = [1e0,1e0,1e0,1e0]
var_scale_label = [L"\mathbf{\times 10^{8}}",L"\mathbf{\times 10^{2}}",L"\mathbf{\times 10^{2}}",L"\mathbf{\times 10}",] 

# Loop over the coefficients
for n in 1:n_coeff
        # Setup the ticks in the axes
        mean_axs[n] = set_ticks(mean_axs[n], μ[1:(end-1)], coeff_mean[:,n]./mean_scale_value[n], n_ticks = 1)
        var_axs[n] = set_ticks(var_axs[n], μ[1:(end-1)], coeff_var[:,n]./var_scale_value[n], n_ticks = 1)

        # Plot the data
        lines!(mean_axs[n], μ[1:(end-1)], coeff_mean[:,n]./mean_scale_value[n], linewidth = 5, color = (:blue,0.5))
        lines!(var_axs[n], μ[1:(end-1)], coeff_var[:,n]./var_scale_value[n], linewidth = 5, color = (:red,0.5))

        # Add the labels to show the scale for the y-axes
        #Label(fig[n, 1, Top()], halign = :left, mean_scale_label[n], fontsize = 50)
        #Label(fig[n, 2, Top()], halign = :right, var_scale_label[n], fontsize = 50)
end

# Export the figure
save("../fig/coefficients.png", fig)
