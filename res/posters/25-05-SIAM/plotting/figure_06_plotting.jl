include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
equilibria = readin("../data/figure_06/equilibria.csv")
μ = equilibria[:,1]
eq = equilibria[:,2]

# Import a reppresentative solution in the ensemble
sol = readin("../data/figure_06/solutions/1.csv")
Ne = length(sol[:,1]) 

# Import the Taylor coefficients
taylor_coefficients = readin("../data/figure_06/taylor.csv")
Nμ = length(taylor_coefficients[:,1])

# Import the solutions of the ensemble mean non-linear least-squares problem
c_nonlinear = readin("../data/figure_06/mean_coefficients.csv")
Nc = length(c_nonlinear[1,:])

# Import the mean guess for the coefficients
c_linear = readin("../data/figure_06/mean_guess.csv") 

# Import the mean escape rates
escapes = readin("../data/figure_06/escape_rates.csv") 
escape_analytic = escapes[:,1]
escape_linear = escapes[:,2]
escape_nonlinear = escapes[:,3]

# Import the escape ratios
escape_ratios = readin("../data/figure_06/escaped_percentage.csv")

# Import error estimates
estimates = readin("../data/figure_06/error_estimates.csv") 
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

# Array of indices to export the plot
idx = [1, Nμ]

# Define counting index
global cnt = 1

#############################
# Coefficient C3 statistics #
#############################
        
# Create and customise the C3 coefficient statistics 
global fig, AX = mkfig(size = [1200,2400],
                box_position = [3:4,1:2],
                limits = ((μ[1],μ[end]), (0.8,1.25)),
                lab = [L"\mathbf{\mu}", L"c^{*}_3"],
                lab_size = [60,70],
                lab_pad = [-60.0,-90.0],
                x_ticks = [μ[1],μ[end]],
                y_ticks = [0.75,1.50],
                ticks_lab_trunc = [1,1]
)
# Plot the heatmap associated to the distribution of the non-linear solutions
for n in 1:Nμ
        #heatmap!(ax, ones(Nb).*μ[n], bins[:,3,n], pdf[:,3,n], colormap = Reverse(:deep))
end

# Loop over the parameter values
@showprogress for n in idx 
        # Get the parameter values
        μ0 = trunc(μ[n], digits=1)
        # Define boolean variables for the plot formats
        y_lab = true 
        y_ticks_lab = true 
        x_ticks_lab_pos = :right
 
        # Assign the values of the boolean variables according to the counting index
        if cnt == 2
                y_lab = false 
                y_ticks_lab = false 
                x_ticks_lab_pos = :left
        end
        
        # Check if least one particle survived (didn't escape the basin)
        if escape_ratios[n] < 1.0::Float64
                # Get the current parameter index as valid
                push!(valid, n) 

                # Import the indices the unescaped trajectories
                unescaped = readin("../data/figure_06/escapes/$n.csv")
      
                #########################
                # Approximate histogram #
                #########################

                # Define the plot limits for the equilibrium distribution histogram 
                x_inf = -(1.75::Float64)
                x_sup = 4.75::Float64
                y_inf = -(0.25::Float64)
                y_sup = 6.5::Float64 

                # Define the domain of the non-linear least-squares solution for the stationary distribution
                x_uns = -(1/(3*c_nonlinear[n,3]))*(sqrt((c_nonlinear[n,2])^2 - 3*c_nonlinear[n,1]*c_nonlinear[n,3]) + c_nonlinear[n,2])
                domain = LinRange(x_uns, x_sup, 1000) 

                # Create and customise the histogram figure 
                global fig, ax = mkfig(fig = fig,
                                title = L"\text{\mu=%$(μ0)}",
                                toggle_title = true,
                                title_size = 50,
                                title_gap = 14.0,
                                box_position = [1,cnt],
                                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                                toggle_lab = [false,y_lab],
                                lab = [L"\mathbf{x}", L"\mathbf{p_s(\bar{x}_t)}"],
                                lab_size = [60,60],
                                lab_pad = [-60.0,-70.0],
                                x_ticks = [x_inf,x_sup],
                                y_ticks = [0,y_sup],
                                toggle_ticks_lab = [false,y_ticks_lab],
                                ticks_lab_trunc = [0,1]
                               )
                # Loop over the ensemble particles
                for m in unescaped 
                        # Import the data
                        distribution = readin("../data/figure_06/distribution/$n/$m.csv")
                        # Plot the histogram approximating the stationary distribution 
                        lines!(ax, distribution[:,1], distribution[:,2], color = (:darkgoldenrod2,0.05), linewidth = 6)
                end
                # Overlay the solution of the ensemble mean of the non-linear least-squares problem
                lines!(ax, domain, [p(x, c_nonlinear[n,:]) for x in LinRange(x_uns, x_sup, 1000)], linewidth = 5, color = (:teal,1.0))
                lines!(ax, [x_inf,x_uns], [0,0], linewidth = 5, color = (:teal,1.0))

                ####################
                # Scalar potential #
                ####################
 
                # Define the analytic potential and its linear and non-linear least-squares approximation 
                U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4
                V_linear(x) = c_linear[n,1]*x + c_linear[n,2]*(x^2) + c_linear[n,3]*(x^3)
                V_nonlinear(x) = c_nonlinear[n,1]*x + c_nonlinear[n,2]*(x^2) + c_nonlinear[n,3]*(x^3)

                # Define the limits for the scalar potential plot
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
                global fig, ax = mkfig(fig = fig,
                                box_position = [2,cnt],
                                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                                toggle_lab = [true,y_lab],
                                lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                                lab_pad = [-60.0,-80.0],
                                x_ticks = [x_inf,x_sup],
                                y_ticks = [y_inf,y_sup],
                                toggle_ticks_lab = [true,y_ticks_lab],
                                ticks_lab_xpos = [x_ticks_lab_pos,:top],
                                ticks_lab_trunc = [0,1]
                               )
                # Import the non-linear least-squares coefficients of the entire ensemble
                coefficients = readin("../data/figure_06/fit/$n.csv")
                # Plot the scalar potential
                lines!(ax, domain, [U(x) for x in domain], color = (:darkgoldenrod2,1.00), linewidth = 8)
                lines!(ax, domain, [V_nonlinear(x) for x in domain] .- shift, color = (:teal,1.00), linewidth = 5)
        end

        # Update the plotting counter
        global cnt = cnt + 1
end

# Plot the real value of c3 
lines!(AX, μ, taylor_coefficients[:,1], color = (:darkgoldenrod2,1.00), linewidth = 8)
# Plot the ensemble mean non-linear solution for the C3 
#lines!(AX, μ[valid], c_nonlinear[valid,3], color = (:teal,1.00), linewidth = 5)
lines!(AX, μ, c_nonlinear[:,3], color = (:teal,1.00), linewidth = 5)

# Export the coefficients statistics plot 
save("../fig/fig6.png", fig)
