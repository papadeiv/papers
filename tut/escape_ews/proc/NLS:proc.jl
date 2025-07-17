include("../../../inc/IO.jl")
include("../../../inc/EscapeProblem.jl")
include("../../../inc/SystemAnalysis.jl")
include("../../../inc/PotentialLearning.jl")
include("../../../inc/TimeseriesAnalysis.jl")
include("../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Define the noise level
σ = 0.250::Float64

# Define the standard deviation of the perturbation of the initial guess
β = 1e-6

# Define the width of the sliding window and the split point 
width = 0.300::Float64
check = 0.100::Float64

# Find the tipping point
solution = readin("../data/NLS/solution.csv")
tip_chk, tip_idx = find_tipping(solution[:,3], width, check)
writeout([tip_idx, t[tip_idx], μ[tip_idx]], "../data/NLS/tipping_point.csv")

# Extract the subseries prior the tipping point
t = solution[1:tip_idx,1]
μ = solution[1:tip_idx,2]
u = solution[1:tip_idx,3]
Nt = length(t)

# Assemble the sliding window
window = get_window_parameters(Nt, width)
Nw = window[1]
Ns = window[2]

# Compute the variance EWSs
μ_var, u_var = ews.variance(μ, u, width)

# Number of bins for the histograms across the sliding window (2% of the number of timesteps in the window)
Nb = convert(Int64, floor(0.02*Nw))

# Define the number of coefficients for the non-linear solution
Nc = convert(Int64, 3e0)

# Define the deterministic dynamics (drift) of the fast variable
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define empty matrix to store the solutions of the optimisation problem
coefficients = Matrix{Float64}(undef, Ns, Nc)

# Define empty arrays to store the escape rate ews
u_esc = Vector{Float64}(undef, Ns)
true_esc = Vector{Float64}(undef, Ns)

# Define empty arrays to store the horizontal and vertical shifts 
x_s = Vector{Float64}(undef, Ns)
y_s = Vector{Float64}(undef, Ns)

# Define empty array to store the numerical error 
error = Vector{Float64}(undef, Ns)

# Define empty matrix to store the coefficients of the Taylor's expansion
taylor = Matrix{Float64}(undef, Ns, 4)

# Define empty matrix to store the analysis parameters
parameters = Matrix{Float64}(undef, Ns, 8)

println("Nt = ", Nt, ", Ns = ", Ns)
println("Nw = ", Nw, ", Nb = ", Nb)

# Loop over the window's strides
printstyled("Computing the escape EWS across the sliding window\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Ns
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

        #println("n = ", n, ", n_end = ", n_end)

        # Extract the subseries in the window
        t_t = t[n:n_end]
        μ_t = μ[n:n_end]
        u_t = u[n:n_end]

        # Define empty arrays for the stable and unstable equilibria
        stable = Float64[]
        unstable = Float64[]

        # Loop over the parameter values across the window
        for μt in μ_t
                # Get the equilibria at the current parameter value
                equilibria = get_equilibria(f, μt, domain=[-10,10])
                # Update the array for the stable equilibrium
                push!(stable, sort(equilibria[1], rev=true)[1])
                # Update the array for the unstable equilibrium
                if length(equilibria[2]) > 0
                        push!(unstable, equilibria[2][1])
                end
        end

        # Detrend the timeseries by the drift of the stable equilibrium 
        u_det = detrend(t_t, u_t, alg = "exact", qse = stable)
        trend = u_det[1]
        residuals = u_det[2]
 
        # Fit an empirical distribution to the detrended timeseries
        bins, pdf = fit_distribution(residuals, n_bins=Nb+1)
       
        # Define empty arrays to store the solutions of the LLS and NLLS problems
        guess = Vector{Float64}(undef, Nc)

        # Select an initial guess
        if n == 1
                # Compute the inverted distribution for the scalar potential under OUP assumption
                local xs, Vs = invert_equilibrium_distribution(bins, pdf, σ)

                # Solve the linear least-squares problem to get an initial guess
                guess = (approximate_potential(xs, Vs, degree=Nc).coeffs)[2:(Nc+1)] + β.*randn(Nc)
        else
                # Select the non-linear solution at the previous iteration
                guess = coefficients[(n-1),:] + β.*randn(Nc)
        end

        # Solve the non-linear least-squares problem to fit the scalar potential 
        coefficients[n,:] = fit_potential(bins, pdf, σ, initial_guess=guess)

        # Define the ground truth and its derivatives
        U(x) = μ[n_end]*x + x^2 - x^3 + (1/5)*(x^4)
        Ux(x) = μ[n_end] + 2*x - 3*(x^2) + (4/5)*(x^3)
        Uxx(x) = 2 - 6*x + (12/5)*(x^2) 
        U3x(x) = - 6 + (24/5)*x 

        # Store the coefficients of the Taylor expansion
        taylor[n,1] = U(stable[end]) 
        taylor[n,2] = Ux(stable[end]) 
        taylor[n,3] = Uxx(stable[end])/2
        taylor[n,4] = U3x(stable[end])/6

        # Compute the horizontal shift
        c = coefficients[n,:] 
        bounds = get_stationary_points(Polynomial(prepend!(c, 0.0)))
        x_s[n] = stable[end] - bounds[2]

        # Compute the vertical shift
        y_s[n] = U(stable[end]) - (Polynomial(c))(bounds[2])

        # Define the reconstructed potential and its second derivative
        V(x) = coefficients[n,1]*x + coefficients[n,2]*(x^2) + coefficients[n,3]*(x^3)
        Vxx(x) = 2*coefficients[n,2] + 6*coefficients[n,3]*x

        # Compute the escape rate ews
        true_esc[n] = kramer_escape(U, Uxx, stable[end], unstable[end], σ)
        u_esc[n] = kramer_escape(V, Vxx, bounds[2], bounds[1], σ)

       # Store the analysis parameters
        parameters[n,1] = Uxx(stable[end])
        parameters[n,2] = Uxx(unstable[end])
        parameters[n,3] = U(unstable[end])-U(stable[end])
        parameters[n,4] = exp(-2*((U(unstable[end])-U(stable[end]))/σ^2))
        parameters[n,5] = Vxx(bounds[2])
        parameters[n,6] = Vxx(bounds[1])
        parameters[n,7] = V(bounds[1])-V(bounds[2])
        parameters[n,8] = exp(-2*((V(bounds[1])-V(bounds[2]))/σ^2))

        # Define the shifted potential for computing the error 
        Vs(x) = y_s[n] + coefficients[n,1]*(x - x_s[n]) + coefficients[n,2]*(x - x_s[n])^2 + coefficients[n,3]*(x - x_s[n])^3

        # Define the bounds of integration for the error
        true_coefficients = [0.0, μ[n_end], +1.0, -1.0, 0.2]
        bounds = get_stationary_points(Polynomial(true_coefficients))
        domain = LinRange(bounds[2], bounds[3], 1000) 

        # Compute the L-2 norm of the reconstruction error between the ground truth and the numerical approximation 
        error[n] = sqrt(sum([(U(x) - Vs(x))^2 for x in domain])) 

        # Export the data 
        writeout(hcat(t_t, μ_t, trend), "../data/NLS/trend/$n.csv")
        writeout(residuals, "../data/NLS/residuals/$n.csv")
        writeout(hcat(bins, pdf), "../data/NLS/histograms/$n.csv")
end

# Export the solutions of the non-linear least-squares problem
writeout(coefficients, "../data/NLS/coefficients.csv")

# Export the EWSs
writeout(hcat(μ_var, u_var, u_esc, true_esc), "../data/NLS/ews.csv")

# Export the horizontal and vertical shifts
writeout(hcat(x_s, y_s), "../data/NLS/shifts.csv")

# Export the numerical error
writeout(error, "../data/NLS/error.csv")

# Export the Taylor's expansion
writeout(taylor, "../data/NLS/taylor.csv")

# Export the analysis parameters 
writeout(parameters, "../data/NLS/parameters.csv")
