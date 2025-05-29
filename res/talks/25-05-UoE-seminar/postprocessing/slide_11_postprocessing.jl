include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Define the noise-level of the SDE
σ = 0.250::Float64

# Define the standard deviation of the perturbation of the initial guess
η = 1e-2

# Import the data from csv
solution = readin("../data/slides_10_to_13/solution.csv")
t = solution[:,1]
Nt = length(t)

residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])
Ns = (Nt - Nw) + 1

# Dimension of the parameter space of the optimisation problem 
Nc = convert(Int64,3e0)

# Matrices to store the solutions {c0, c1, c2, c3} of the linear and non-linear least-squares problems
mean_linear_coefficients = Matrix{Float64}(undef, Ns, Nc)
mean_nonlinear_coefficients = Matrix{Float64}(undef, Ns, Nc)
lin_linear_coefficients = Matrix{Float64}(undef, Ns, Nc)
lin_nonlinear_coefficients = Matrix{Float64}(undef, Ns, Nc)

# Loop over the window's strides
printstyled("Finding the local minima of the non-linear optimisation problem\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Ns
        # Import the windowed histogram 
        mean_hist = readin("../data/slide_10/histograms/SN_mean/$n.csv")
        mean_bins = mean_hist[:,1]
        mean_pdf = mean_hist[:,2]
        lin_hist = readin("../data/slide_10/histograms/SN_lin/$n.csv")
        lin_bins = lin_hist[:,1]
        lin_pdf = lin_hist[:,2]

        ############################
        # Mean detrended histogram #
        ############################
        
        # Select an initial guess
        initial_guess = Vector{Float64}(undef, Nc)
        if n == 1
                # Compute the inverted distribution for the scalar potential under OUP assumption
                local xs, Vs = invert_equilibrium_distribution(mean_bins, mean_pdf, σ)

                # Solve the linear least-squares problem to get an initial guess
                linear_solution = approximate_potential(xs, Vs, degree=Nc).coeffs
                initial_guess = linear_solution[2:end] + η.*randn(Nc)
                mean_linear_coefficients[n,:] = initial_guess
        else
                # Select the non-linear solution at the previous iteration
                initial_guess = mean_nonlinear_coefficients[(n-1),:] + η.*randn(Nc)
        end

        # Solve the non-linear least-squares problem to fit the scalar potential 
        optimal_solution = fit_potential(mean_bins, mean_pdf, σ, initial_guess=initial_guess)
        mean_nonlinear_coefficients[n,:] = optimal_solution 

        ##############################
        # Linear detrended histogram #
        ##############################
        
        # Select an initial guess
        if n == 1
                # Compute the inverted distribution for the scalar potential under OUP assumption
                local xs, Vs = invert_equilibrium_distribution(lin_bins, lin_pdf, σ)

                # Solve the linear least-squares problem to get an initial guess
                linear_solution = approximate_potential(xs, Vs, degree=Nc).coeffs
                initial_guess = linear_solution[2:end] + η.*randn(Nc)
                lin_linear_coefficients[n,:] = initial_guess
        else
                # Select the non-linear solution at the previous iteration
                initial_guess = lin_nonlinear_coefficients[(n-1),:] + η.*randn(Nc)
        end

        # Solve the non-linear least-squares problem to fit the scalar potential 
        optimal_solution = fit_potential(lin_bins, lin_pdf, σ, initial_guess=initial_guess)
        lin_nonlinear_coefficients[n,:] = optimal_solution
end

# Export the solutions of the optimisation problem 
writeout(mean_linear_coefficients, "../data/slide_11/mean_linear_coefficients.csv")
writeout(mean_nonlinear_coefficients, "../data/slide_11/mean_nonlinear_coefficients.csv")
writeout(lin_linear_coefficients, "../data/slide_11/lin_linear_coefficients.csv")
writeout(lin_nonlinear_coefficients, "../data/slide_11/lin_nonlinear_coefficients.csv")
