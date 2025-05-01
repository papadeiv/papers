include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
equilibria = readin("../data/equilibria.csv")
μ = equilibria[:,1]
eq = equilibria[:,2]

# Import a reppresentative solution in the ensemble
sol = readin("../data/solutions/1.csv")
Ne = length(sol[:,1]) 

# Import the solutions of the ensemble mean non-linear least-squares problem
c_nonlinear = readin("../data/mean_coefficients.csv")
Nc = length(c_nonlinear[1,:])
Nμ = length(c_nonlinear[:,1])

# Import the mean guess for the coefficients
c_linear = readin("../data/mean_guess.csv") 

# Import the mean escape rates
escapes = readin("../data/escape_rates.csv") 
escape_analytic = escapes[:,1]
escape_linear = escapes[:,2]
escape_nonlinear = escapes[:,3]

# Import the escape ratios
escape_ratios = readin("../data/escaped_percentage.csv")

# Import error estimates
estimates = readin("../data/error_estimates.csv") 
error_linear = estimates[:,1]
error_nonlinear = estimates[:,2]

# Number of bins for the histogram of ensemble coefficients 
Nb = convert(Int64,3e1)

# Define the stochastic diffusion
σ = 0.200::Float64
D = (σ^2)/2.0::Float64

# Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
c0(μ) = - μ[1]*xs(μ) - μ[2]*(xs(μ))^2 - μ[3]*(xs(μ))^3

# Define an arbitrary cubic with the the above constraint on {c0}
V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
 
# Define the exact stationary probability distribution
f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))
N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf), parameters=μ)
p(x, μ) = N(μ)*f(x, μ)

# Array to store the indices of the parameter values for which at least one trajectory in the ensemble survived
valid = Int64[]

# Loop over the parameter values
printstyled("Generating figures of the non-linear optimisation\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Check if least one particle survived (didn't escape the basin)
        if escape_ratios[n] < 1.0::Float64
                # Get the current parameter index as valid
                push!(valid, n) 

                # Import the indices the unescaped trajectories
                unescaped = readin("../data/escapes/$n.csv")
      
                #########################
                # Approximate histogram #
                #########################

                # Define the plot limits for the equilibrium distribution histogram 
                x_inf = -(2.00::Float64)
                x_sup = 4.00::Float64
                y_inf = -(0.25::Float64)
                y_sup = 6.25::Float64 

                # Define the domain of the non-linear least-squares solution for the stationary distribution
                x_uns = -(1/(3*c_nonlinear[n,3]))*(sqrt((c_nonlinear[n,2])^2 - 3*c_nonlinear[n,1]*c_nonlinear[n,3]) + c_nonlinear[n,2])
                domain = LinRange(x_uns, x_sup, 1000) 

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
                # Loop over the ensemble particles
                for m in unescaped 
                        # Import the data
                        distribution = readin("../data/distribution/$n/$m.csv")
                        # Plot the histogram approximating the stationary distribution 
                        lines!(ax, distribution[:,1], distribution[:,2], color = (:black,0.15), linewidth = 1)
                end
                # Overlay the solution of the ensemble mean of the non-linear least-squares problem
                lines!(ax, domain, [p(x, c_nonlinear[n,:]) for x in LinRange(x_uns, x_sup, 1000)], linewidth = 4, color = (:blue,0.3))

                # Export the equilibrium distribution histogram
                save("../fig/distribution/$n.png", fig)

                ####################
                # Scalar potential #
                ####################
 
                # Define the analytic potential and its linear and non-linear least-squares approximation 
                U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4
                V_linear(x) = c_linear[n,1]*x + c_linear[n,2]*(x^2) + c_linear[n,3]*(x^3)
                V_nonlinear(x) = c_nonlinear[n,1]*x + c_nonlinear[n,2]*(x^2) + c_nonlinear[n,3]*(x^3)

                # Define the limits for the scalar potential plot
                x_inf = -(1.75::Float64)
                x_sup = 4.75::Float64
                y_inf = -(1.0::Float64)
                y_sup = 3.5::Float64 

                # Define the domain of the function
                domain = LinRange(x_inf, x_sup, 1000)

                # Find optimal shift for the linear-least squares guess 
                x_stb_g = (1/(3*c_linear[n,3]))*(sqrt((c_linear[n,2])^2 - 3*c_linear[n,1]*c_linear[n,3]) - c_linear[n,2])
                shift_g = abs(V_linear(x_stb_g)-U(eq[n]))

                # Find optimal shift for the non linear-least squares solution 
                x_stb = (1/(3*c_nonlinear[n,3]))*(sqrt((c_nonlinear[n,2])^2 - 3*c_nonlinear[n,1]*c_nonlinear[n,3]) - c_nonlinear[n,2])
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
                # Import the non-linear least-squares coefficients of the entire ensemble
                coefficients = readin("../data/fit/$n.csv")
                # Plot the scalar potential
                lines!(ax, domain, [U(x) for x in domain], color = (:black,0.25), linewidth = 4.5)
                # Plot the mean polynomial least-squares solution
                lines!(ax, domain, [V_linear(x) for x in domain] .- shift_g, color = (:red,0.50), linewidth = 3)
                lines!(ax, domain, [V_nonlinear(x) for x in domain] .- shift, color = (:blue,0.50), linewidth = 3)
 
                # Export the scalar potential function plot 
                save("../fig/fit/$n.png", fig)
        end
end

#############################
# Ensemble mean escape rate #
#############################

# Create and customise the ensemble mean escape rate figure 
fig, ax = mkfig(size = [1800,1800],
                box_position = [1,1],
                bg_out = "#eeeeeeff",
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\textbf{mean escape rate}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                ax_scale = [identity, identity],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0,0.15],
                toggle_ticks_lab = [false,true],
                ticks_lab_trunc = [0,2]
)
# Plot the analytic escape rate 
lines!(ax, μ[valid], escape_analytic[valid], color = (:black,0.25), linewidth = 10)
# Plot the escape rate for the linear reconstruction 
lines!(ax, μ[valid], escape_linear[valid], color = (:red,1.00), linewidth = 6)
# Plot the escape rate for the non-linear reconstruction 
lines!(ax, μ[valid], escape_nonlinear[valid], color = (:blue,1.00), linewidth = 3)

# Add the escape percentage axis 
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\textbf{% escapes}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0,1],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,1]
)
scatter!(ax, μ[valid], escape_ratios[valid], color = (:darkgreen,0.50), markersize = 40, strokecolor = (:black,1.00), strokewidth = 3)
 
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
# Plot the heatmap associated to the distribution of the non-linear solutions
for n in 1:Nμ
        #heatmap!(ax, ones(Nb).*μ[n], bins[:,1,n], pdf[:,1,n], colormap = Reverse(:deep))
end
# Plot the initial guess for the current coefficient 
lines!(ax, μ[valid], c_linear[valid,1], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ[valid], c_nonlinear[valid,1], color = (:blue,1.00), linewidth = 4)
 
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
# Plot the heatmap associated to the distribution of the non-linear solutions
for n in 1:Nμ
        #heatmap!(ax, ones(Nb).*μ[n], bins[:,2,n], pdf[:,2,n], colormap = Reverse(:deep))
end
# Plot the initial guess for the current coefficient 
lines!(ax, μ[valid], c_linear[valid,2], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ[valid], c_nonlinear[valid,2], color = (:blue,1.00), linewidth = 4)
 
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
# Plot the heatmap associated to the distribution of the non-linear solutions
for n in 1:Nμ
        #heatmap!(ax, ones(Nb).*μ[n], bins[:,3,n], pdf[:,3,n], colormap = Reverse(:deep))
end
# Plot the initial guess for the current coefficient 
lines!(ax, μ[valid], c_linear[valid,3], color = (:red,1.00), linewidth = 4)
# Plot the ensemble mean non-linear solution for the current coefficient
lines!(ax, μ[valid], c_nonlinear[valid,3], color = (:blue,1.00), linewidth = 4)

# Export the coefficients statistics plot 
save("../fig/ensemble_error_statistics.png", fig)

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
lines!(ax, μ[valid], error_linear[valid], color = (:red,1.00), linewidth = 4)
lines!(ax, μ[valid], error_nonlinear[valid], color = (:blue,1.00), linewidth = 4)

# Add the escape percentage axis 
fig, ax = mkfig(fig=fig,
                box_position = [2,1],
                limits = ((μ[1],μ[end]), nothing),
                lab = [L"\mathbf{\mu}", L"\textbf{% escapes}"],
                toggle_lab = [true,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0,1],
                toggle_ticks_lab = [true,true],
                ticks_lab_trunc = [1,1]
)
scatter!(ax, μ[valid], escape_ratios[valid], color = (:darkgreen,0.50), markersize = 40, strokecolor = (:black,1.00), strokewidth = 3)

# Export the error analysis plot 
save("../fig/ensemble_error_analysis.png", fig)
