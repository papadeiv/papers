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

# Number of observations per sample 
N_observations = convert(Int64,1e1)

# Standard deviation of the observations 
σ = 0.100

# Generate normally distributed random samples
ξ = σ.*randn(N_observations)
        
# Generate the observations
x = LinRange(domain[1], domain[end], N_observations)
y = [f(xn) for xn in x] # Deterministic data
z = y + ξ               # Noisy data

# Fit the noisy data and get its residual 
U, r = fit_potential(x, z, degree=d)
# Fit the deterministic data (residual will be 0)
V = fit_potential(x, y, degree=d)[1]

# Create the figure for th samples, the analytcial function and the least-square solution
fig, ax = mkfig()
set_ticks(ax, x, y, n_ticks = 3)

# Plot the fit of the deterministic data
lines!(ax, LinRange(domain[1], domain[end], 1000), [V(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 5, color = (:red, 1.0))
# Plot the fit of the noisy data 
lines!(ax, LinRange(domain[1], domain[end], 1000), [U(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = (:blue, 0.5), linestyle = :dash)
# Plot the analytical function
lines!(ax, LinRange(domain[1], domain[end], 1000), [f(x) for x in LinRange(domain[1], domain[end], 1000)], linewidth = 3, color = (:black, 0.5))
# Plot the vertical lines from the samples to the least-square solution using the residual
for m in 1:N_observations
        lines!(ax, [x[m], x[m]], [z[m], z[m]-r[m]], linewidth = 1, color = :black)
end
# Plot the samples
scatter!(ax, x, z, markersize = 20, color = :blue, strokewidth = 3, strokecolor = :black)

# Export the figure
save("../fig/quadratic_fit.png", fig)
