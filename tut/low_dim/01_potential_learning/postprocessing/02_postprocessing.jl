include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the parameters values
μ = readin("../data/parameter.csv")
# Import the solutions
data = readin("../data/solution.csv") 
t = data[:,1]
u = data[:,2:end]

# Define the exact scalar potential function 
U(x, μ) = μ*x + x^2 - x^3 + (1/5)*x^4
# Define the noise-level of the SDE
σ = 0.200::Float64

# Number of fixed values
Nμ = length(μ)
# Number of timesteps
Nt = length(t)
# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e2)
# Number of polynomial coefficients of the least-squares regression
Nc = convert(Int64,4e0)

# Empty array to store the normalisation constants
N = Vector{Float64}(undef, Nμ)
N_OUP = Vector{Float64}(undef, Nμ)

# Equilibria
eq = [2.73045::Float64, 2.72759::Float64]

# Loop over the paramater's values
printstyled("Fitting the data with polynomial least squares\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Fit a normalised empirical distribution to the data
        u_μ = u[:,n]
        bins, pdf = fit_distribution(u_μ, n_bins=Nb+1)

        # Export the normalised distribution data
        writeout(hcat(bins, pdf), "../data/distribution/$n.csv")

        # Approximate the normalisation constant
        p(x, null_p) = exp(-(2*U(x, μ[n]))/(σ^2))
        N[n] = get_normalisation_constant(p, (-Inf, Inf), accuracy=1e-10)

        # Invert the normalised equilibrium distribution using the true normalisation 
        xs, Vs = invert_equilibrium_distribution(bins, pdf, σ, N=N[n])

        # Export the potential data
        writeout(hcat(xs, Vs), "../data/potential/$n.csv") 

        # Invert the normalised equilibrium distribution using the OUP normalisation 
        ϴ = abs(-2 + 6*eq[n] - (12/5)*(eq[n])^2)
        N_OUP[n] = sqrt(ϴ/(pi*(σ^2))) 
        xs, Vs = invert_equilibrium_distribution(bins, pdf, σ, N=N_OUP[n])
        writeout(hcat(xs, Vs), "../data/potential/OUP$n.csv") 
       
        # Approximate the potential function using the empirical distribution
        V = approximate_potential(xs, Vs, degree=Nc-1)
        coefficients = V.coeffs

        # Export the coefficients data
        writeout(coefficients, "../data/fit/$n.csv")
end

# Export the normalisation constants
writeout(hcat(N, N_OUP), "../data/normalisation_constants.csv")
