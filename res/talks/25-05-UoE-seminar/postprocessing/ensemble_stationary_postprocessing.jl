include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the parameters values
equilibria = readin("../data/equilibria.csv")
μ = equilibria[:,1]
Nμ = length(μ)

# Import a reppresentative solution in the ensemble
sol = readin("../data/solutions/1.csv")
Ne = length(sol[:,1]) 

# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e2)
# Dimension of the parameter space of the optimisation problem 
Nc = convert(Int64,3e0)

# Define the dynamics
f(x, μ) = -μ - 2*x + 3*(x^2) -(4/5)*x^3
# Define the noise-level of the SDE
σ = 0.200::Float64

# Matrix to store the ensemble mean linear (guess) and non-linear least-squares solutions 
mean_guess = Matrix{Float64}(undef, Nμ, Nc)
mean_coefficients = Matrix{Float64}(undef, Nμ, Nc)

# Vectors to store the error estimates
error_linear = Vector{Float64}(undef, Nμ)
error_nonlinear = Vector{Float64}(undef, Nμ)

# Vectors to store the mean escape rates 
escape_analytic = Vector{Float64}(undef, Nμ)
escape_linear = Vector{Float64}(undef, Nμ)
escape_nonlinear = Vector{Float64}(undef, Nμ)

# Vector to store the percentage of particles of the ensemble that escape from the basin
escaped = Vector{Float64}(undef, Nμ)

# Loop over the paramater's values
printstyled("Finding the local minima of the non-linear optimisation problem\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Create the directory to store the histograms
        mkpath("../data/distribution/$n/")
        # Import the ensemble solution at the current parameter value μ
        u_μ = readin("../data/solutions/$n.csv")

        # Compute the unstable equilibria at the current parameter value
        xu = sort(get_equilibria(f, μ[n], domain=[-20,20])[2], rev=true)[1]

        # Array to store the outcomes of the check for escapes for the particles in the ensemble
        unescaped = Int64[] 
        # Matrix to store the solutions {c0, c1, c2, c3} of the linear least-squares
        initial_guesses = Matrix{Float64}(undef, Ne, Nc+1)

        # Loop over the ensemble solutions (optimal initial guess for the non-linear regression)
        for m in 1:Ne
                # Extract the trajectory of the m-th particle
                local u = u_μ[m,:]

                # Check if the trajectory has escaped the basin of attraction
                check = check_escapes(u, [xu, Inf])

                # Proceed if the trajectory has not escaped
                if !check
                        # Fit an empirical distribution to the data
                        bins, pdf = fit_distribution(u, n_bins=Nb+1)
                        # Export the histogram of the current particle 
                        writeout(hcat(bins, pdf), "../data/distribution/$n/$m.csv")

                        # Compute the inverted distribution for the scalar potential under OUP assumption
                        xs, Vs = invert_equilibrium_distribution(bins, pdf, σ)
                        # Derive the coefficients via linear least-squares
                        initial_guesses[m,:] = approximate_potential(xs, Vs, degree=3).coeffs

                        # Update the check index array
                        push!(unescaped, m)
                end
        end

        # Export the escape indices array 
        writeout(unescaped, "../data/escapes/$n.csv")

        # Compute the number of trajectories that haven't escaped
        N_unescaped = length(unescaped)
        escaped[n] = (Ne - N_unescaped)/Ne

        # Proceed if at least one trajectory has remained bounded
        if N_unescaped > 0
                # Compute the ensemble mean for the initial guess without using {c0}
                mean_guess[n,:] = [mean(initial_guesses[[unescaped[k] for k in 1:N_unescaped], m]) for m in 2:(Nc+1)]

                # Matrix to store the solutions {c1, c2, c3} of the non-linear least squares regression 
                coefficients = Matrix{Float64}(undef, N_unescaped, Nc)

                # Local index for storing the results
                local idx = 1
 
                # Loop over the unescaped ensemble solutions (actual non-linear regression)
                for m in unescaped 
                        # Extract the trajectory of the m-th particle
                        local u = u_μ[m,:]
                        # Fit an empirical distribution to the data
                        bins, pdf = fit_distribution(u, n_bins=Nb+1)
                        # Non-linear least-squares fit of the coefficients of the scalar potential 
                        coefficients[idx,:] = fit_potential(bins, pdf, σ, initial_guess=mean_guess[n,:])
                        # Update the local index
                        idx = idx + 1
                end
                # Export the coefficients matrix
                writeout(coefficients, "../data/fit/$n.csv")

                # Compute the ensemble mean coefficients
                mean_coefficients[n,:] = [mean(coefficients[:,m]) for m in 1:Nc]

                # Matrices for the ensemble distribution of coefficients
                coefficients_bins = Matrix{Float64}(undef, Nb, Nc)
                coefficients_pdf = Matrix{Float64}(undef, Nb, Nc)

                # Fit an histogram to the ensemble distribution of coefficients
                for m in 1:Nc
                        coefficients_bins[:,m], coefficients_pdf[:,m] = fit_distribution(coefficients[:,m], n_bins=Nb+1)
                end

                # Export the ensemble distribution of coefficients
                writeout(coefficients_bins, "../data/coefficients/bins$n.csv")
                writeout(coefficients_pdf, "../data/coefficients/pdf$n.csv")

                # Define the analytic potential and its linear and non-linear least-squares approximations...
                U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4
                V_linear(x) = mean_guess[n,1]*x + mean_guess[n,2]*(x^2) + mean_guess[n,3]*(x^3)
                V_nonlinear(x) = mean_coefficients[n,1]*x + mean_coefficients[n,2]*(x^2) + mean_coefficients[n,3]*(x^3)
                # ... and their second derivatives
                Uxx(x) = 2.0::Float64 - 6.0::Float64 + (12.0::Float64/5.0::Float64)*x^2
                Vxx_linear(x) = 2.0::Float64*mean_guess[n,2] + 6.0::Float64*mean_guess[n,3]*x
                Vxx_nonlinear(x) = 2.0::Float64*mean_coefficients[n,2] + 6.0::Float64*mean_coefficients[n,3]*x

                # Define the bounds of integration for the error
                true_coefficients = [0.0, μ[n], +1.0, -1.0, 0.2]
                bounds = get_stationary_points(Polynomial(true_coefficients))
                domain = LinRange(bounds[2], bounds[3], 1000) 

                # Compute the L-2 norm of the error of the linear and non-linear least-squares solutions w.r.t. the analytic potential 
                error_linear[n] = sqrt(sum([(U(x) - V_linear(x))^2 for x in domain])) 
                error_nonlinear[n] = sqrt(sum([(U(x) - V_nonlinear(x))^2 for x in domain])) 

                # Define the local extrema of the 3 potentials
                a = bounds[2] 
                b = bounds[3]
                a_linear = -(1/(3*mean_guess[n,3]))*(sqrt((mean_guess[n,2])^2 - 3*mean_guess[n,1]*mean_guess[n,3]) + mean_guess[n,2]) 
                b_linear = +(1/(3*mean_guess[n,3]))*(sqrt((mean_guess[n,2])^2 - 3*mean_guess[n,1]*mean_guess[n,3]) - mean_guess[n,2])
                a_nonlinear = -(1/(3*mean_coefficients[n,3]))*(sqrt((mean_coefficients[n,2])^2 - 3*mean_coefficients[n,1]*mean_coefficients[n,3]) + mean_coefficients[n,2])
                b_nonlinear = +(1/(3*mean_coefficients[n,3]))*(sqrt((mean_coefficients[n,2])^2 - 3*mean_coefficients[n,1]*mean_coefficients[n,3]) - mean_coefficients[n,2])

                # Compute the mean escape rates of the 3 potentials 
                escape_analytic[n] = kramer_escape(U, Uxx, a, b, σ)
                escape_linear[n] = kramer_escape(V_linear, Vxx_linear, a_linear, b_linear, σ)
                escape_nonlinear[n] = kramer_escape(V_nonlinear, Vxx_nonlinear, a_nonlinear, b_nonlinear, σ)
        end
end

# Export the ensemble mean linear (guess) and non-linear least-squares solutions
writeout(mean_guess, "../data/mean_guess.csv")
writeout(mean_coefficients, "../data/mean_coefficients.csv")

# Export the error estimates
writeout(hcat(error_linear, error_nonlinear), "../data/error_estimates.csv")

# Export the escape rates
writeout(hcat(escape_analytic, escape_linear, escape_nonlinear), "../data/escape_rates.csv")

# Export the percentage of escapes
writeout(escaped, "../data/escaped_percentage.csv")

#=
u_μ = readin("../data/solutions/9.csv")
xu = sort(get_equilibria(f, μ[9], domain=[-20,20])[2], rev=true)[1]
include("../../../../inc/PlottingTools.jl")
fig, ax = mkfig(size = [1200,800],
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
)
lines!(ax, [1,length(u_μ[1,:])], [xu,xu], color = (:blue, 1.00), linewidth = 3, linestyle = :dash)
for n in 1:Ne
        if u_μ[n,end] < xu 
                lines!(ax, LinRange(1,length(u_μ[n,:]),length(u_μ[n,:])), u_μ[n,:], color = (:brown2, 0.35), linewidth = 1)
        else
                lines!(ax, LinRange(1,length(u_μ[n,:]),length(u_μ[n,:])), u_μ[n,:], color = (:black, 0.15), linewidth = 1)
        end
end
save("../fig/escapes.png", fig)
=#
