include("../../../inc/IO.jl")
include("../../../inc/EscapeProblem.jl")
include("../../../inc/SystemAnalysis.jl")
include("../../../inc/PotentialLearning.jl")
include("../../../inc/TimeseriesAnalysis.jl")

# Define the noise level and the bifurcation parameter
σ = 0.125::Float64
D = (σ^2)/2
μ = 1.200::Float64

# Define the standard deviation of the perturbation of the initial guess
β = 1e-3

# Define the ground truth and its derivatives
U(x) = μ*x + x^2 - x^3 + (1/5)*(x^4)
Ux(x) = μ + 2*x - 3*(x^2) + (4/5)*(x^3)
Uxx(x) = 2 - 6*x + (12/5)*(x^2) 
U3x(x) = - 6 + (24/5)*x 

# Define the stable and unstable equilibria of the ground truth
u_stable = 2.52::Float64
u_unstable = 2.52::Float64

# Compute and export the coefficients of the Taylor expansion of the ground truth around the stable equilibria
ct1 = Ux(u_stable)
ct2 = Uxx(u_stable)/2
ct3 = U3x(u_stable)/6
writeout(hcat(ct1, ct2, ct3), "../data/stationary_ensemble/taylor.csv")

# Import the data from csv 
solutions = readin("../data/stationary_ensemble/solutions.csv")
t = readin("../data/stationary_ensemble/time.csv")
Ne = length(solutions[:,1])
Nt = length(solutions[1,:])

# Number of bins for the histograms (2% of the number of timesteps of each solution)
Nb = convert(Int64, floor(0.02*Nt))

# Define the number of coefficients for the non-linear solution
Nc = convert(Int64, 3e0)

# Define empty matrix to store the solutions of the optimisation problem
coefficients = Matrix{Float64}(undef, Ne, Nc)

# Define empty matrix for the statistical analysis of the coefficients
parameters = Matrix{Float64}(undef, Ne, 4)

# Define empty arrays to store the escape rate ews
u_esc = Vector{Float64}(undef, Ne)
true_esc = Vector{Float64}(undef, Ne)

# Define empty arrays to store the horizontal and vertical shifts 
x_s = Vector{Float64}(undef, Ne)
y_s = Vector{Float64}(undef, Ne)

# Define empty array to store the numerical error 
error = Vector{Float64}(undef, Ne)

println("Nb = ", Nb)

# Loop over the window's strides
printstyled("Performing statistical analysis on the least-squares solutions\n"; bold=true, underline=true, color=:light_blue)
#=@showprogress=#for n in 1:Ne
        # Extract the current solution from the ensemble and center it
        u = solutions[n,:] .- u_stable 

        # Fit an empirical distribution to the timeseries
        bins, pdf = fit_distribution(u, n_bins=Nb+1)
       
        # Compute the inverted distribution under OUP assumption
        xs, Vs = invert_equilibrium_distribution(bins, pdf, σ)

        # Solve the linear least-squares problem to get an initial guess
        guess = (approximate_potential(xs, Vs, degree=Nc).coeffs)[2:(Nc+1)] + β.*randn(Nc)
        #guess = [ct1, ct2, ct3]
        guess = [0.0, 1.0, 0.1] + β.*randn(Nc)
        #guess = β.*randn(Nc)
        println(n, ")  ", guess)

        # Solve the non-linear least-squares problem to fit the potential 
        coefficients[n,:] = fit_potential(bins, pdf, σ, initial_guess=guess)
        c1 = coefficients[n,1] 
        c2 = coefficients[n,2] 
        c3 = coefficients[n,3] 

        # Define the reconstructed potential and its second derivative
        V(x) = c1*x + c2*(x^2) + c3*(x^3)
        Vxx(x) = 2*c2 + 6*c3*x

        # Compute the random variables targeted by the analysis 
        a = (1/(3*c3))*(sqrt((c2)^2 - 3*c1*c3) - c2)      
        parameters[n,1] = a               # Estimated stable equilibrium
        parameters[n,2] = V(a)            # Potential value at equilibrium
        parameters[n,3] = Vxx(a)          # Curvature value at equilibrium
        parameters[n,4] = exp(V(a)/D)     # Escape probability at equilibrium

        # Compute the horizontal shift
        c = coefficients[n,:] 
        bounds = get_stationary_points(Polynomial(prepend!(c, 0.0)))
        x_s[n] = u_stable - bounds[2]

        # Compute the vertical shift
        y_s[n] = U(u_stable) - (Polynomial(c))(bounds[2])

        # Compute the escape rate ews
        true_esc[n] = kramer_escape(U, Uxx, u_stable, u_unstable, σ)
        u_esc[n] = kramer_escape(V, Vxx, bounds[2], bounds[1], σ)

        # Define the shifted potential for computing the error 
        Vs(x) = y_s[n] + c1*(x - x_s[n]) + c2*(x - x_s[n])^2 + c3*(x - x_s[n])^3

        # Define the bounds of integration for the error
        true_coefficients = [0.0, μ, +1.0, -1.0, 0.2]
        bounds = get_stationary_points(Polynomial(true_coefficients))
        domain = LinRange(bounds[2], bounds[3], 1000) 

        # Compute the L-2 norm of the reconstruction error between the ground truth and the numerical approximation 
        error[n] = sqrt(sum([(U(x) - Vs(x))^2 for x in domain])) 
end

# Export the solutions of the non-linear least-squares problem
writeout(coefficients, "../data/stationary_ensemble/coefficients.csv")

# Export the EWSs
writeout(hcat(u_esc, true_esc), "../data/stationary_ensemble/ews.csv")

# Export the horizontal and vertical shifts
writeout(hcat(x_s, y_s), "../data/stationary_ensemble/shifts.csv")

# Export the numerical error
writeout(error, "../data/stationary_ensemble/error.csv")

# Export the analysis parameters 
writeout(parameters, "../data/stationary_ensemble/parameters.csv")
