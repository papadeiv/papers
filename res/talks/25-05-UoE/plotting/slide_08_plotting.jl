include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level
σ = 0.250::Float64
D = (σ^2)/(2.0::Float64) 

# Import the data 
escapes = readin("../data/slide_08/escapes.csv")
N = escapes[:,1]
tip_prob = escapes[:,2]
solution = readin("../data/slide_08/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nμ = length(μ)

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)
# Loop over the realisations
@showprogress for n in 1:Nμ
        # Extract and truncate the tipping probability
        tip = trunc(100*tip_prob[n], digits = 5)

        ######################
        # Potential function #
        ######################
 
        # Define the analytic potential function
        V(x) = μ[n]*x + x^2 - x^3 + (1/5)*(x^4)

        # Define the domain for plotting the potential
        x_inf = -(2.0::Float64)
        x_sup = 4.0::Float64
        domain = LinRange(x_inf,x_sup,1000)

        # Create and customise the potential axis 
        fig, ax = mkfig(size = [1800,1200],
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [2,1:2],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((x_inf, x_sup), (-2, 4)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        lab_pad = [-60.0,-50.0],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [-2,4],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the scalar potential 
        lines!(ax, domain, [V(x) for x in domain], linewidth = 4.5, color = (:black, 1.0))
        # Plot the solution at the current timestep
        scatter!(ax, u[n], V(u[n]), color = (:teal, 1.0), markersize = 30, strokewidth = 3)

        ############################
        # Equilibrium distribution #
        ############################
        
        # Define the normalised equilibrium density
        ρ(x) = N[n]*exp(-(1.0::Float64/D)*(V(x)))

        # Create and customise the density axis 
        fig, ax = mkfig(fig = fig,
                        #bg_out = colorant"rgb(147,237,213)",
                        box_position = [1,1:2],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((x_inf, x_sup), (-0.25, 6)),
                        lab = [L"\mathbf{x}", L"\mathbf{\rho(x)}"],
                        toggle_lab = [false,true],
                        lab_pad = [-60.0,-30.0],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [0,6],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        toggle_ticks_lab = [false,true],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the equilibrium distribution 
        lines!(ax, domain, [ρ(x) for x in domain], linewidth = 4.5, color = (:black, 1.0))
        # Plot a thick line showing the valued of the density at the current timestep
        lines!(ax, [u[n],u[n]], [0, ρ(u[n])], linewidth = 8, color = (:teal, 0.25))
        # Plot the solution at the current timestep
        scatter!(ax, u[n], 0, color = (:teal, 1.0), markersize = 30, strokewidth = 3)
        # Print the value of the tipping probability at the current timestep
        text!(u[n], 5.7, text = L"\mathbf{p}_{\textbf{tip}}=\mathbf{%$tip%}", color = (:teal, 1.00), align = (:center,:center), fontsize = 30)

        #################################
        # Escape probability timeseries #
        #################################

        # Create and customise the escape axis 
        fig, ax = mkfig(fig = fig,
                        box_position = [1,3:5],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((t[1], t[end]), (-0.05, 1.05)),
                        lab = [L"\textbf{time}", L"\mathbf{p}_{\textbf{tip}}"],
                        toggle_lab = [false,true],
                        lab_pad = [-60.0,-30.0],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        x_ticks = [t[1],t[end]],
                        y_ticks = [0,1],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        toggle_ticks_lab = [false,true],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the escape probability timeseries 
        lines!(ax, t[1:n], tip_prob[1:n], linewidth = 4.5, color = (:teal, 1.0))
        # Plot the solution at the current timestep
        scatter!(ax, t[n], tip_prob[n], color = (:teal, 1.0), markersize = 30, strokewidth = 3)
 
        #######################
        # Solution timeseries #
        #######################

        # Create and customise the solution axis 
        fig, ax = mkfig(fig = fig,
                        box_position = [2,3:5],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((t[1], t[end]), (x_inf, x_sup)),
                        lab = [L"\textbf{time}", L"\mathbf{x}"],
                        lab_pad = [-60.0,-50.0],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        x_ticks = [t[1],t[end]],
                        y_ticks = [x_inf,x_sup],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the escape probability timeseries 
        lines!(ax, t[1:n], u[1:n], linewidth = 2.5, color = (:teal, 1.0))
        # Plot the solution at the current timestep
        scatter!(ax, t[n], u[n], color = (:teal, 1.0), markersize = 30, strokewidth = 3)
 
        # Export the figure
        save("../fig/slide_08/$n.png", fig)
end
