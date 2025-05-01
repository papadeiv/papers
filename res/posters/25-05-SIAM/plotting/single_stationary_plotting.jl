include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
equilibria = readin("../data/equilibria.csv")
μ = equilibria[:,1]
eq = equilibria[:,2]

# Import the linear and nonlinear solutions for the coefficients 
nonlinear = readin("../data/nonlinear.csv")
linear = readin("../data/linear.csv") 
Nc = length(nonlinear[1,:])
Nμ = length(nonlinear[:,1])

# Import the mean escape rates
escapes = readin("../data/escape_rates.csv") 
escape_analytic = escapes[:,1]
escape_linear = escapes[:,2]
escape_nonlinear = escapes[:,3]

# Import the escape ratios
escape_time = readin("../data/escape_time.csv")

# Import the variance EWS 
variance = readin("../data/variance.csv")

# Import error estimates
estimates = readin("../data/error_estimates.csv") 
error_linear = estimates[:,1]
error_nonlinear = estimates[:,2]

# Empty arrays to store the error in the pdfs
pdf_error_linear = Vector{Float64}(undef, Nμ)
pdf_error_nonlinear = Vector{Float64}(undef, Nμ)

# Define the stochastic diffusion
σ = 0.200::Float64
D = (σ^2)/2.0::Float64

# Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
c0(μ) = - μ[1]*xs(μ) - μ[2]*(xs(μ))^2 - μ[3]*(xs(μ))^3

# Define an arbitrary cubic with the the above constraint on {c0}
V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
 
# Define the stationary probability distribution
f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))
N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf), parameters=μ)
p(x, μ) = N(μ)*f(x, μ)

# Loop over the parameter values
printstyled("Generating figures of the non-linear optimisation\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        #########################
        # Approximate histogram #
        #########################

        # Define the plot limits for the equilibrium distribution histogram 
        x_inf = -(2.00::Float64)
        x_sup = 4.00::Float64
        y_inf = -(0.25::Float64)
        y_sup = 6.25::Float64 

        # Define the domain of the non-linear least-squares solution for the stationary distribution
        x_uns = -(1/(3*nonlinear[n,3]))*(sqrt((nonlinear[n,2])^2 - 3*nonlinear[n,1]*nonlinear[n,3]) + nonlinear[n,2])
        domain = LinRange(x_uns, x_sup, 1000) 

        # Define the analytic pdf and its linear and non-linear reconstructions
        p_linear(x) = p(x, linear[n,:])
        p_nonlinear(x) = p(x, nonlinear[n,:])

        # Create and customise the histogram figure 
        fig, ax = mkfig(size = [1200,900],
                                bg_out = "#eeeeeeff",
                                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                                lab = [L"\mathbf{x}", L"\textbf{density}"],
                                lab_pad = [-60.0,-40.0],
                                x_ticks = [x_inf,x_uns,x_sup],
                                y_ticks = [0,y_sup],
                                ticks_lab_trunc = [1,1]
                               )
        # Import the data
        distribution = readin("../data/distribution/$n.csv")
        # Plot the histogram approximating the stationary distribution 
        lines!(ax, distribution[:,1], distribution[:,2], color = (:black,0.15), linewidth = 1)
        # Overlay the solution of the linear least-squares problem
        lines!(ax, domain, [p_linear(x) for x in domain], linewidth = 3, color = (:brown2,1.00))
        # Overlay the solution of the non-linear least-squares problem
        lines!(ax, domain, [p_nonlinear(x) for x in domain], linewidth = 1.5, color = (:blue,0.75))

        # Export the equilibrium distribution histogram
        save("../fig/distribution/$n.png", fig)

        # Compute the L2-norm of the error of the pdfs
        pdf_error_linear[n] = norm([(p_linear(distribution[m,1]) - distribution[m,2]) for m in 1:length(distribution[:,1])] ,2)
        pdf_error_nonlinear[n] = norm([(p_nonlinear(distribution[m,1]) - distribution[m,2]) for m in 1:length(distribution[:,1])] ,2)

        ####################
        # Scalar potential #
        ####################
 
        # Define the analytic potential and its linear and non-linear least-squares approximation 
        U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4
        V_linear(x) = linear[n,1]*x + linear[n,2]*(x^2) + linear[n,3]*(x^3)
        V_nonlinear(x) = nonlinear[n,1]*x + nonlinear[n,2]*(x^2) + nonlinear[n,3]*(x^3)

        # Define the limits for the scalar potential plot
        x_inf = -(1.75::Float64)
        x_sup = 4.75::Float64
        y_inf = -(1.0::Float64)
        y_sup = 3.5::Float64 

        # Define the domain of the function
        domain = LinRange(x_inf, x_sup, 1000)

        # Find optimal shift for the linear-least squares guess 
        x_stb_g = (1/(3*linear[n,3]))*(sqrt((linear[n,2])^2 - 3*linear[n,1]*linear[n,3]) - linear[n,2])
        shift_g = abs(V_linear(x_stb_g)-U(eq[n]))

        # Find optimal shift for the non linear-least squares solution 
        x_stb = (1/(3*nonlinear[n,3]))*(sqrt((nonlinear[n,2])^2 - 3*nonlinear[n,1]*nonlinear[n,3]) - nonlinear[n,2])
        shift = abs(V_nonlinear(x_stb)-U(eq[n]))

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
        lines!(ax, domain, [U(x) for x in domain], color = (:black,0.25), linewidth = 4.5)
        # Plot the mean polynomial least-squares solution
        lines!(ax, domain, [V_linear(x) for x in domain] .- shift_g, color = (:red,0.50), linewidth = 3)
        lines!(ax, domain, [V_nonlinear(x) for x in domain] .- shift, color = (:blue,0.50), linewidth = 3)
 
        # Export the scalar potential function plot 
        save("../fig/fit/$n.png", fig)
end

###############
# Escape rate #
###############

# Create and customise the escape rate figure 
fig, ax = mkfig(size = [1800,1200],
                bg_out = "#eeeeeeff",
                box_position = [1,1],
                limits = ((μ[1],μ[end]), (-0.01, 0.15)),
                lab = [L"\mathbf{\mu}", L"\textbf{escape rate}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                ax_scale = [identity, identity],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0,0.15],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,2]
)
# Plot the analytic escape rate 
lines!(ax, μ, escape_analytic, color = (:black,1.00), linewidth = 6)
# Plot the escape rate for the linear reconstruction 
lines!(ax, μ, escape_linear, color = (:red,0.75), linewidth = 3)
# Plot the escape rate for the non-linear reconstruction 
lines!(ax, μ, escape_nonlinear, color = (:blue,0.75), linewidth = 3)
# Plot the variance of the trajectory at each parameter value
lines!(ax, μ, variance, color = (:lime,1.00), linewidth = 3)

#=
# Create a mirrored axis for the variance
mirror_ax = mirror_axis(fig, [μ[1],μ[end]], y_lab = L"\textbf{variance}", color = :darkgreen)
# Plot the variance of the trajectory at each parameter value
lines!(mirror_ax, μ, variance, color = (:lime,1.00), linewidth = 3)
=#

# Add the escape time index axis 
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\textbf{escape time}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0, escape_time[1]],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,0]
)
scatter!(ax, μ, escape_time, color = (:darkgreen,0.50), markersize = 40, strokecolor = (:black,1.00), strokewidth = 3)
 
# Export the escape rate figure
save("../fig/escape_rate.png", fig)

#############################
# Coefficient C1 statistics #
#############################
        
# Create and customise the ensemble error statistics figure 
fig, ax = mkfig(size = [1200,1800],
                bg_out = "#eeeeeeff",
                box_position = [1,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\mathbf{\bar{c}_1}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [6,14],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [1,1]
)
# Customise ticks
set_ticks(ax, μ, linear[:,1])
# Plot the initial guess for the current coefficient 
lines!(ax, μ, linear[:,1], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ, nonlinear[:,1], color = (:blue,1.00), linewidth = 4)
 
#############################
# Coefficient C2 statistics #
#############################
        
# Create and customise the ensemble error statistics figure 
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\mathbf{\bar{c}_2}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [-7,-3],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [1,1]
)
# Customise ticks
set_ticks(ax, μ, linear[:,2])
# Plot the initial guess for the current coefficient 
lines!(ax, μ, linear[:,2], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ, nonlinear[:,2], color = (:blue,1.00), linewidth = 4)
 
#############################
# Coefficient C3 statistics #
#############################
        
# Create and customise the ensemble error statistics figure 
fig, ax = mkfig(fig=fig,
                box_position = [3,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\mathbf{\bar{c}_3}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0.5,1.00],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,2]
)
# Customise ticks
set_ticks(ax, μ, linear[:,3])
# Plot the initial guess for the current coefficient 
lines!(ax, μ, linear[:,3], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ, nonlinear[:,3], color = (:blue,1.00), linewidth = 4)

# Export the coefficients statistics plot 
save("../fig/coefficients.png", fig)

##################
# Error analysis #
##################
 
# Create and customise the error decay figure
fig, ax = mkfig(size = [1800,1800],
                box_position = [1,1],
                limits = ((μ[1],μ[end]), nothing),
                bg_out = "#eeeeeeff",
                lab = [L"\mathbf{\mu}", L"\mathbf{||V-U||_2}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [50,300],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [1,0]
)
# Customise ticks
set_ticks(ax, μ, error_linear)
# Plot the reconstruction error from linear regression
lines!(ax, μ, error_linear, color = (:red,1.00), linewidth = 4)
# Plot the reconstruction error from non-linear regression
lines!(ax, μ, error_nonlinear, color = (:blue,1.00), linewidth = 4)

# Add the escape percentage axis 
#=
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\textbf{escape time}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0, escape_time[1]],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,0]
)
scatter!(ax, μ, escape_time, color = (:darkgreen,0.50), markersize = 40, strokecolor = (:black,1.00), strokewidth = 3)
=#
# Add the error in the pdf axis
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\mathbf{||p-q||_2}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                #y_ticks = [0, escape_time[1]],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,2]
)
# Customise ticks
set_ticks(ax, μ, pdf_error_linear)
# Plot the reconstruction error from linear regression
lines!(ax, μ, pdf_error_linear, color = (:red,1.00), linewidth = 4)
# Plot the reconstruction error from non-linear regression
lines!(ax, μ, pdf_error_nonlinear, color = (:blue,1.00), linewidth = 4)

# Export the error analysis plot 
save("../fig/enhanced_error_analysis.png", fig)
