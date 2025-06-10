include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

printstyled("Postprocessing the data\n"; bold=true, underline=true, color=:darkgreen)

# Define the fixed parameter value
μ0 = 0.625::Float64

# Define the noise level
σ = 0.125::Float64

# Read in the timeseries of the stationary process
solution = readin("../data/slide_05/solution.csv")
t = solution[:,1]
u = solution[:,2]
Nt = length(t)

# Define number of bins for the histogram (0.25% of the total number of timesteps)
Nb = convert(Int64, floor(0.00125*Nt))

# Define the number of coefficients for the non-linear solution
Nc = convert(Int64, 3e0)

# Build a histogram of the process
bins, pdf = fit_distribution(u, n_bins=Nb+1)

# Find an optimal initial guess for the non-linear solver
xs, Vs = invert_equilibrium_distribution(bins, pdf, σ)
guess = (approximate_potential(xs, Vs, degree=Nc).coeffs)[2:(Nc+1)]

# Compute the non-linear solution
coefficients = fit_potential(bins, pdf, σ, initial_guess=guess)

# Export the data so far
writeout(hcat(bins, pdf), "../data/slide_05/histogram.csv")
writeout(hcat(guess, coefficients), "../data/slide_05/coefficients.csv")

# Compute a vertical shift to match the true potential and the reconstruction 
true_coefficients = [0.0, μ0, +1.0, -1.0, 0.2]
true_bounds = get_stationary_points(Polynomial(true_coefficients))
U(x) = (Polynomial(true_coefficients))(x)
V(x) = (Polynomial(coefficients))(x)
bounds = get_stationary_points(Polynomial(prepend!(coefficients, 0.0)))
ys = U(true_bounds[3]) - V(bounds[2])

# Export the shift 
writeout(hcat(ys, ys), "../data/slide_05/shift.csv")
