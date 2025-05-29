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

# Import the variance EWS 
stats = readin("../data/statistics.csv")

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
printstyled("Generating the figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Nμ
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
