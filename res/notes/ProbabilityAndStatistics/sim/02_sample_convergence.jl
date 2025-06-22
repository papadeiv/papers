include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")

# Define the domain of the observations
a = -1.000 
b = -a 
domain = [a, b]

# Define the deterministic function
f(x) = x^2

# define the degree of the polynomial least-squares regression
d = 2

# Vector of different number of observations per sample 
N_observations = convert(Int64,1e1)

# Standard deviation of the observations 
σ = 0.100

# Number of samples
N_samples = convert(Int64,1e1)

# Create the figure for the mean coefficients reconstruction 
global Fig, Ax = mkfig()

# Loop over the number of samples
using ProgressMeter
@showprogress for j in 1:N_samples
        # Generate random samples
        z = σ.*randn(N_observations)
        
        # Generate noisy function 
        x = LinRange(domain[1], domain[end], N_observations)
        y = [f(xn) for xn in x] + z 

        # Fit a parabola using the implemented least-squares fit
        U, C = fit_potential(x, y, degree=d)

        # Extract the diagonal of the covariance matrix
        covs = [C[n] for n in findall(x -> x > 0.0, Diagonal(C))]

        # Store the least-squares coefficients of the fit
        coefficients[j,:] = U.coeffs

        # Create the figure for th samples, the analytcial function and the least-square solution
        local fig, ax = mkfig()
        set_ticks(ax, x, y, n_ticks = 3)

        # Plot the implemented least-squares fit
        lines!(Ax, LinRange(domain[1], domain[end], 1000), [U(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = (:black, 0.1))
end

# Compute the mean and std of each coefficient
mean_coeff = [mean(coefficients[:,n]) for n in 1:(d+1)]
std_coeff = [std(coefficients[:,n]) for n in 1:(d+1)]

# Construct the mean polynomial
U_m(x) = Polynomial(mean_coeff)(x) 

# Construct the std polynomial
U_s_pos(x) = Polynomial(mean_coeff + std_coeff)(x) 
U_s_neg(x) = Polynomial(mean_coeff - std_coeff)(x) 

# Plot the analytical function 
lines!(Ax, LinRange(domain[1], domain[end], 1000), [f(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 5, color = :red)
# Plot the mean least-squares fit
lines!(Ax, LinRange(domain[1], domain[end], 1000), [U_m(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 5, color = (:blue, 1.0), linestyle = :dash)
# Plot the std least-squares fit
lines!(Ax, LinRange(domain[1], domain[end], 1000), [U_s_pos(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = (:green, 1.0))
lines!(Ax, LinRange(domain[1], domain[end], 1000), [U_s_neg(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = (:green, 1.0))

# Export the figure
save("../fig/mean_fit.png", Fig)
