include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Define the noise level
σ = 0.125::Float64

# Define the standard deviation of the perturbation of the initial guess
β = 1e-5

# Import the data from csv 
solution = readin("../data/slide_06/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(t)

# Define the width of the sliding window (as a relative length of the timeseries)
width = 0.200::Float64

# Assemble the sliding window
window = get_window_parameters(Nt, width)
Nw = window[1]
Ns = window[2]

# Compute the variance EWSs
t_var, u_var = ews.variance(t, u, width)

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

display(Nt)
display(Ns)
display(Nw)
display(Nb)

# Loop over the window's strides
printstyled("Computing the escape EWS across the sliding window\n"; bold=true, underline=true, color=:light_blue)
#=@showprogress=#for n in 1:Ns
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

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

        # Define a boolean and an integer value for tipping detection
        tipped = false
        tip_idx = 0::Int64
        
        # Check if you are in a region of bistability
        if length(unstable) > 0
                # Check the presence of a tipping point
                tipped, tip_idx = check_escapes(u_t, [unstable[end], +Inf])
        end

        # Define empty arrays for the trend and the residuals in the window 
        trend = Vector{Float64}(undef, Nw)
        residuals = Vector{Float64}(undef, Nw)

        # Check if the tipping point has been reached
        if tipped
                # Export the tipping point and stop the loop
                writeout([n, (n + tip_idx - 1), t_t[tip_idx], μ_t[tip_idx]], "../data/slide_06/tipping_point.csv")
                break
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

        # Define the true potential function and its second derivative
        U(x) = μ[n_end]*x + x^2 - x^3 + (1/5)*x^4
        Uxx(x) = 2 - 6*x +(12/5)*(x^2) 

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

        display(n)
        display(μ_t[end])
        display(u_esc[n])
        display(true_esc[n])

        # Export the data 
        writeout(hcat(t_t, μ_t, trend), "../data/slide_06/trend/$n.csv")
        writeout(residuals, "../data/slide_06/residuals/$n.csv")
        writeout(hcat(bins, pdf), "../data/slide_06/histograms/$n.csv")
end

# Export the solutions of the non-linear least-squares problem
writeout(coefficients, "../data/slide_06/coefficients.csv")

# Export the EWSs
writeout(hcat(t_var, u_var, u_esc, true_esc), "../data/slide_06/ews.csv")

# Export the horizontal and vertical shifts
writeout(hcat(x_s, y_s), "../data/slide_06/shifts.csv")
