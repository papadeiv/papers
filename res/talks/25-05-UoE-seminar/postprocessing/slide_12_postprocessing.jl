include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Import the data from csv
mean_nonlinear_coefficients = readin("../data/slide_11/mean_nonlinear_coefficients.csv")
lin_nonlinear_coefficients = readin("../data/slide_11/lin_nonlinear_coefficients.csv")
Ns = length(lin_nonlinear_coefficients[:,1])
Nc = length(lin_nonlinear_coefficients[1,:])

solution = readin("../data/slides_10_to_13/solution.csv")
μ = solution[:,2]
Nt = length(μ)

residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])

# Arrays to store the shifts 
xs_mean = Vector{Float64}(undef, Ns)
xs_lin = Vector{Float64}(undef, Ns)
ys_mean = Vector{Float64}(undef, Ns)
ys_lin = Vector{Float64}(undef, Ns)

# Arrays to store the error estimates 
x_uns = Vector{Float64}(undef, Ns)
x_stb = Vector{Float64}(undef, Ns)
mean_error = Vector{Float64}(undef, Ns)
lin_error = Vector{Float64}(undef, Ns)

# Loop over the window's strides
printstyled("Computing the L2-norm of the reconstruction error\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Ns
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

        # Extract the current solutions
        c_mean = prepend!(mean_nonlinear_coefficients[n,:], 0.0)
        c_lin = prepend!(lin_nonlinear_coefficients[n,:], 0.0)

        # Define the analytic potential function
        U(x) = μ[n_end]*x + x^2 - x^3 + (1/5)*x^4

        # Define the bounds of integration for the error
        true_coefficients = [0.0, μ[n_end], +1.0, -1.0, 0.2]
        bounds = get_stationary_points(Polynomial(true_coefficients))
        x_uns[n] = bounds[2]
        x_stb[n] = bounds[3]
        domain = LinRange(x_uns[2], x_stb[3], 5000) 

        # Compute the horizontal shift
        mean_bounds = get_stationary_points(Polynomial(c_mean))
        xs_mean[n] = x_stb[n] - mean_bounds[2]
        lin_bounds = get_stationary_points(Polynomial(c_lin))
        xs_lin[n] = x_stb[n] - lin_bounds[2]

        # Compute the vertical shift
        V_mean(x) = (Polynomial(c_mean))(x)
        ys_mean[n] = U(x_stb[n]) - V_mean(mean_bounds[2])
        V_lin(x) = (Polynomial(c_lin))(x)
        ys_lin[n] = U(x_stb[n]) - V_lin(lin_bounds[2])

        # Define the new, shifted reconstructions of the potential
        Vs_mean(x) = ys_mean[n] + c_mean[2]*(x-xs_mean[n]) + c_mean[3]*(x-xs_mean[n])^2 + c_mean[4]*(x-xs_mean[n])^3
        Vs_lin(x) = ys_lin[n] + c_lin[2]*(x-xs_lin[n]) + c_lin[3]*(x-xs_lin[n])^2 + c_lin[4]*(x-xs_lin[n])^3

        # Compute the L-2 norm of the reconstruction error
        mean_error[n] = norm([U(x) - Vs_mean(x) for x in domain], 2)
        lin_error[n] = norm([U(x) - Vs_lin(x) for x in domain], 2)
end

# Export the shifts 
writeout(hcat(xs_mean, xs_lin, ys_mean, ys_lin), "../data/slide_12/shifts.csv")
# Export the error estimates 
writeout(hcat(μ[Nw:end], x_uns, x_stb, mean_error, lin_error), "../data/slide_12/errors.csv")
