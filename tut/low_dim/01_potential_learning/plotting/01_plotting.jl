include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/IO.jl")

# Define the exact scalar potential function 
U(x, μ) = μ*x + (1/3)*x^3

# Import the parameter range
μ = readin("../data/parameter.csv")
# Get the number of parameter's values
Nμ = length(μ)

# Loop over the parameter values
printstyled("Creating the plots\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        ####################### 
        # Stationary solution #
        ####################### 
 
        # Import the data
        solution = readin("../data/solutions/monostable_$n.csv")
        histogram = readin("../data/distribution/monostable_$n.csv")
        coefficients = readin("../data/coefficients/monostable_$n.csv")
        potential = readin("../data/potential/monostable_$n.csv")

        # Define the polynomial least-squares solution
        V = Polynomial(coefficients)

        # Define the limits for the potential plot
        x_inf = -(3.0::Float64)
        x_sup = 3.0::Float64
        y_inf = -(4.0::Float64)
        y_sup = 7.0::Float64 
        # Define the domain of the function
        domain = LinRange(x_inf, x_sup, 1000)
        # Create the figure for the potential reconstruction
        fig, ax = mkfig(size = [1000, 1000],
                        bg_out = :white,
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        toggle_lab = [true, true],
                        lab_pad = [-60.0, -60.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [y_inf,y_sup],
                        toggle_ticks_lab = [true, true],
                        ticks_lab_trunc = [1,1]
                       )
        # Plot the scalar potential 
        lines!(ax, domain, [U(x, μ[n]) for x in domain], color = (:black,0.25), linewidth = 4.5)
        # Plot the mean polynomial least-squares solution
        lines!(ax, domain, [V(x) for x in domain], color = (:red,1.00), linewidth = 4)
        # Plot the datapoints from the inverted distribution
        scatter!(ax, potential[:,1], potential[:,2], marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Export the figure
        save("../fig/potential/$n.png", fig)

        # Create the figure for the full solution timeseries 
        fig, ax = mkfig(size = [1600, 1600],
                        bg_out = :white,
                        pad = (100,60,10,20),
                        box_position = [1,1:3],
                        lab = [L"\mathbf{t}", L"\mathbf{u}"],
                        toggle_lab = [false, true],
                        lab_pad = [-60.0, -150.0],
                        toggle_ticks_lab = [false, true],
                        ticks_lab_trunc = [0,4]
                       )
        # Customise the ticks
        set_ticks(ax, solution[:,1], solution[:,2])
        # Plot the timeseries of the full solution
        lines!(ax, solution[:,1], solution[:,2], color = (:red, 0.75), linewidth = 1)

        # Create the figure for the full solution histogram 
        nullfig, ax = mkfig(fig = fig, 
                            box_position = [1,4],
                            lab = [L"\mathbf{p(x)}", L"\mathbf{x}"],
                            toggle_lab = [false, false],
                            lab_pad = [-60.0, -150.0],
                            toggle_ticks_lab = [true, false],
                            ticks_lab_trunc = [1,1]
                           )
        # Customise the ticks
        set_ticks(ax, histogram[:,2], histogram[:,1])
        # Plot the empirical distribution
        scatter!(ax, histogram[:,2], histogram[:,1], markersize = 20, color = (:red, 0.5), strokecolor = :black, strokewidth = 1) 

        ######################
        # Detrended solution #
        ######################
 
        # Import the data
        residuals = readin("../data/solutions/detrended_$n.csv")
        histogram = readin("../data/distribution/detrended_monostable_$n.csv")
 
        # Create the figure for the residual timeseries
        nullfig, ax = mkfig(fig = fig, 
                            box_position = [2,1:3],
                            lab = [L"\mathbf{t}", L"\mathbf{u_{\textbf{res}}:=u-\bar{u}}"],
                            toggle_lab = [true, true],
                            lab_pad = [-60.0, -150.0],
                            toggle_ticks_lab = [true, true],
                            ticks_lab_trunc = [0,4]
                           )
        # Customise the ticks
        set_ticks(ax, solution[:,1], residuals)
        # Plot the timeseries of the residual 
        lines!(ax, solution[:,1], residuals, color = (:blue, 0.75), linewidth = 1)

        # Create the figure for the residual histogram 
        nullfig, ax = mkfig(fig = fig, 
                            box_position = [2,4],
                            lab = [L"\mathbf{p(x)}", L"\mathbf{x}"],
                            toggle_lab = [false, false],
                            lab_pad = [-60.0, -150.0],
                            toggle_ticks_lab = [true, false],
                            ticks_lab_trunc = [1,1]
                           )
        # Customise the ticks
        set_ticks(ax, histogram[:,2], histogram[:,1])
        # Plot the empirical distribution
        scatter!(ax, histogram[:,2], histogram[:,1], markersize = 20, color = (:blue, 0.5), strokecolor = :black, strokewidth = 1) 

        ################################
        # Topologically equivalent OUP #
        ################################

        # Import the data
        OUP = readin("../data/solutions/OUP_$n.csv")
        histogram = readin("../data/distribution/equivalent_OUP_$n.csv")

        # Create the figure for the topologically equivalent OUP 
        nullfig, ax = mkfig(fig = fig, 
                            box_position = [3,1:3],
                            lab = [L"\mathbf{t}", L"\mathbf{u_{\textbf{OUP}}}"],
                            toggle_lab = [true, true],
                            lab_pad = [-60.0, -150.0],
                            toggle_ticks_lab = [true, true],
                            ticks_lab_trunc = [0,4]
                           )
        # Customise the ticks
        set_ticks(ax, OUP[:,1], OUP[:,2])
        # Plot the timeseries of the OUP 
        lines!(ax, OUP[:,1], OUP[:,2], color = (:green, 0.75), linewidth = 1)

        # Create the figure for the OUP histogram 
        nullfig, ax = mkfig(fig = fig, 
                            box_position = [3,4],
                            lab = [L"\mathbf{p(x)}", L"\mathbf{x}"],
                            toggle_lab = [true, false],
                            lab_pad = [-60.0, -150.0],
                            toggle_ticks_lab = [true, false],
                            ticks_lab_trunc = [1,1]
                           )
        # Customise the ticks
        set_ticks(ax, histogram[:,2], histogram[:,1])
        # Plot the empirical distribution
        scatter!(ax, histogram[:,2], histogram[:,1], markersize = 20, color = (:green, 0.5), strokecolor = :black, strokewidth = 1) 

        # Export the figure
        save("../fig/solutions/$n.png", fig)
end
