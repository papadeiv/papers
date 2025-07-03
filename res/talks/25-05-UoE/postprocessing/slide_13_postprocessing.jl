include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Import the data from csv
mean_nonlinear_coefficients = readin("../data/slide_11/mean_nonlinear_coefficients.csv")
lin_nonlinear_coefficients = readin("../data/slide_11/lin_nonlinear_coefficients.csv")
Ns = length(lin_nonlinear_coefficients[:,1])

solution = readin("../data/slides_10_to_13/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(μ)

residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])

# Compute and export the variance ews
width = 0.150::Float64
μ_ews, var_ews = ews.variance(μ, u, width)
writeout(hcat(μ_ews, var_ews), "../data/slide_13/variance.csv")

# Define the noise-level of the SDE
σ = 0.250::Float64

# Arrays to store the escape rates 
escape_true = Vector{Float64}(undef, Ns)
escape_mean = Vector{Float64}(undef, Ns)
escape_lin = Vector{Float64}(undef, Ns)

# Loop over the window's strides
printstyled("Computing the escape rates\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Ns
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

        # Extract the current solutions
        c_mean = mean_nonlinear_coefficients[n,:]
        c_lin = lin_nonlinear_coefficients[n,:]

        # Define the various potential functions...
        U(x) = μ[n_end]*x + x^2 - x^3 + (1/5)*x^4
        V_mean(x) = c_mean[1]*(x) + c_mean[2]*(x^2) + c_mean[3]*(x^3)
        V_lin(x) = c_lin[1]*(x) + c_lin[2]*(x^2) + c_lin[3]*(x^3)
        #... and their second derivatives
        Uxx(x) = 2 - 6*x + (12/5)*x^2
        Vxx_mean(x) = 2*c_mean[2] + 6*c_mean[3]*x
        Vxx_lin(x) = 2*c_lin[2] + 6*c_lin[3]*x

        # Get the bounds of escape formula
        true_coefficients = [0.0, μ[n_end], +1.0, -1.0, 0.2]
        true_bounds = get_stationary_points(Polynomial(true_coefficients))
        true_b = true_bounds[2]
        true_a = true_bounds[3]
        mean_bounds = get_stationary_points(Polynomial(prepend!(c_mean,0.0)))
        mean_b = mean_bounds[1]
        mean_a = mean_bounds[2]
        lin_bounds = get_stationary_points(Polynomial(prepend!(c_lin,0.0)))
        lin_b = lin_bounds[1]
        lin_a = lin_bounds[2]

        # Compute the escape rates 
        escape_true[n] = kramer_escape(U, Uxx, true_b, true_a, σ)
        escape_mean[n] = kramer_escape(V_mean, Vxx_mean, mean_a, mean_b, σ)
        escape_lin[n] = kramer_escape(V_lin, Vxx_lin, lin_a, lin_b, σ)
end

# Export the escape rates
writeout(hcat(escape_true, escape_mean, escape_lin), "../data/slide_13/escapes.csv")
