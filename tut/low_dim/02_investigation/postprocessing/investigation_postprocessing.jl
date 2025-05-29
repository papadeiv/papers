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
end

# Export the linear (guess) and non-linear least-squares solutions
writeout(linear_coefficients, "../data/linear.csv")
writeout(nonlinear_coefficients, "../data/nonlinear.csv")

# Export the variance
writeout(stats, "../data/statistics.csv")
