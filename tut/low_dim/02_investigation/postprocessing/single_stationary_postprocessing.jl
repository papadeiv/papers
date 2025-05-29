include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the parameters values
equilibria = readin("../data/equilibria.csv")
μ = equilibria[:,1]
Nμ = length(μ)

# Import the solutions
sol = readin("../data/solutions.csv")

# Dimension of the parameter space of the optimisation problem 
Nc = convert(Int64,3e0)

# Define the dynamics
f(x, μ) = -μ - 2*x + 3*(x^2) -(4/5)*x^3
# Define the noise-level of the SDE
σ = 0.200::Float64

# Matrix to store the linear (guess) and non-linear least-squares solutions 
linear_coefficients = Matrix{Float64}(undef, Nμ, Nc)
nonlinear_coefficients = Matrix{Float64}(undef, Nμ, Nc)

# Vectors to store the error estimates
error_linear = Vector{Float64}(undef, Nμ)
error_nonlinear = Vector{Float64}(undef, Nμ)

# Vectors to store the escape rates 
escape_analytic = Vector{Float64}(undef, Nμ)
escape_linear = Vector{Float64}(undef, Nμ)
escape_nonlinear = Vector{Float64}(undef, Nμ)

# Vector to store the escape time index
escape_time = Vector{Int64}(undef, Nμ)

# Vector to store the statistics of the timeseries 
stats = Matrix{Float64}(undef, Nμ, 3)

# Loop over the paramater's values
printstyled("Solving the non-linear optimisation problem using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
#=@showprogress=#Threads.@threads for n in 1:Nμ
        display(μ[n])
        ######################
        # Escape probability #
        ######################
 
        # Import the solution at the current parameter value μ
        u_μ = sol[:,n] 

        # Compute the unstable equilibria at the current parameter value
        xu = sort(get_equilibria(f, μ[n], domain=[-20,20])[2], rev=true)[1]

        # Get the time index of first escape out of the basin of attraction
        check, escape_time[n] = check_escapes(u_μ, [xu, Inf])
        #display(escape_time[n])

        # Define a truncated timeseries of the trajectory until the first escape
        u = u_μ[1:escape_time[n]]

        # Compute the statistics of the truncated timeseries
        stats[n,1] = mean(u)
        stats[n,2] = var(u)
        stats[n,3] = skewness(u)

        # Define number of bins for the histogram (1/50 th of the number of truncated realizations until first exit) to be no less then 5
        Nb = max(5::Int64, convert(Int64, ceil(escape_time[n]/50)))

        # Fit an empirical distribution to the data
        bins, pdf = fit_distribution(u, n_bins=Nb+1)

        # Export the histogram of the current particle 
        writeout(hcat(bins, pdf), "../data/distribution/$n.csv")

        # Compute the inverted distribution for the scalar potential under OUP assumption
        xs, Vs = invert_equilibrium_distribution(bins, pdf, σ)

        # Derive the coefficients via linear least-squares
        linear_coefficients[n,:] = (approximate_potential(xs, Vs, degree=3).coeffs)[2:Nc+1]
        linear_coefficients[n,3] = abs(linear_coefficients[n,3])
        #display(linear_coefficients[n,:])

        # Non-linear least-squares fit of the coefficients of the scalar potential 
        nonlinear_coefficients[n,:] = fit_potential(bins, pdf, σ, initial_guess=linear_coefficients[n,:])

        ##################
        # Error analysis #
        ##################
 
        # Define the analytic potential and its linear and non-linear least-squares approximations...
        U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4
        V_linear(x) = linear_coefficients[n,1]*x + linear_coefficients[n,2]*(x^2) + linear_coefficients[n,3]*(x^3)
        V_nonlinear(x) = nonlinear_coefficients[n,1]*x + nonlinear_coefficients[n,2]*(x^2) + nonlinear_coefficients[n,3]*(x^3)
        # ... and their second derivatives
        Uxx(x) = 2.0::Float64 - 6.0::Float64 + (12.0::Float64/5.0::Float64)*x^2
        Vxx_linear(x) = 2.0::Float64*linear_coefficients[n,2] + 6.0::Float64*linear_coefficients[n,3]*x
        Vxx_nonlinear(x) = 2.0::Float64*nonlinear_coefficients[n,2] + 6.0::Float64*nonlinear_coefficients[n,3]*x

        # Define the bounds of integration for the error
        true_coefficients = [0.0, μ[n], +1.0, -1.0, 0.2]
        bounds = get_stationary_points(Polynomial(true_coefficients))
        domain = LinRange(bounds[2], bounds[3], 1000) 

        # Compute the L-2 norm of the error of the linear and non-linear least-squares solutions w.r.t. the analytic potential 
        error_linear[n] = norm([(U(x) - V_linear(x)) for x in domain], 2)
        error_nonlinear[n] = norm([(U(x) - V_nonlinear(x)) for x in domain], 2)

        # Define the local extrema of the 3 potentials
        a = bounds[2] 
        b = bounds[3]
        a_linear = -(1/(3*linear_coefficients[n,3]))*(sqrt((linear_coefficients[n,2])^2 - 3*linear_coefficients[n,1]*linear_coefficients[n,3]) + linear_coefficients[n,2]) 
        b_linear = +(1/(3*linear_coefficients[n,3]))*(sqrt((linear_coefficients[n,2])^2 - 3*linear_coefficients[n,1]*linear_coefficients[n,3]) - linear_coefficients[n,2])
        a_nonlinear = -(1/(3*nonlinear_coefficients[n,3]))*(sqrt((nonlinear_coefficients[n,2])^2 - 3*nonlinear_coefficients[n,1]*nonlinear_coefficients[n,3]) + nonlinear_coefficients[n,2])
        b_nonlinear = +(1/(3*nonlinear_coefficients[n,3]))*(sqrt((nonlinear_coefficients[n,2])^2 - 3*nonlinear_coefficients[n,1]*nonlinear_coefficients[n,3]) - nonlinear_coefficients[n,2])

        # Compute the mean escape rates of the 3 potentials 
        escape_analytic[n] = kramer_escape(U, Uxx, a, b, σ)
        escape_linear[n] = kramer_escape(V_linear, Vxx_linear, a_linear, b_linear, σ)
        escape_nonlinear[n] = kramer_escape(V_nonlinear, Vxx_nonlinear, a_nonlinear, b_nonlinear, σ)
end

# Export the linear (guess) and non-linear least-squares solutions
writeout(linear_coefficients, "../data/linear.csv")
writeout(nonlinear_coefficients, "../data/nonlinear.csv")

# Export the error estimates
writeout(hcat(error_linear, error_nonlinear), "../data/error_estimates.csv")

# Export the escape rates
writeout(hcat(escape_analytic, escape_linear, escape_nonlinear), "../data/escape_rates.csv")

# Export the escape time index
writeout(escape_time, "../data/escape_time.csv")

# Export the variance
writeout(stats, "../data/statistics.csv")
