include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/IO.jl")

# Number of timesteps of each particle
Nt = convert(Int64,1e5)
# Number of (fixed) parameter values
Nμ = convert(Int64,1e2)
# Number of particles in the ensemble (set this to 1 to reduce the problem to a SDE analysis)
Ne = convert(Int64,1e3)
# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e2)
# Number of polynomial coefficients of the least-squares regression
Nc = convert(Int64,4e0)

# Specify the (additive) noise level
σ = 0.100
g(x) = σ

# Specify the parameter range and export it 
μ = LinRange(-2.000,-0.100,Nμ)
writeout(μ, "../data/parameter.csv")

# Initialise the matrix of polynomial coefficients
coefficients = Matrix{Float64}(undef, Ne, Nc)

# Loop over the parameter range 
printstyled("Simulating the ensemble SDE (monostable saddle node)\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:(Nμ-1)
        # Define the normal form of monostable saddle-node as the first-derivative of its potential 
        U(x, μ) = (1.0/3.0)*x^3 + μ*x
        Ux(x, μ) = x^2 + μ
        Uxx(x, μ) = +2*x

        # Propagate the ensemble forward in time
        t, E = evolve_ensemble(Ux, g, μ[n], Nt=Nt, Ne=Ne)

        # Loop over the particles trajectories
        for m in 1:Ne
                # Extract the trajectory
                u = E[m,:]

                # Get the extrema of one reppresentative trajectory 
                m==1 && global a = minimum(u)
                m==1 && global b = maximum(u)

                # Fit an empirical distribution to the data
                bins, pdf = fit_distribution(u, n_bins=Nb+1)

                # Approximate the potential function using the empirical distribution
                xs, Vs = invert_equilibrium_distribution(bins, pdf, std(u))
                V = approximate_potential(xs, Vs, degree=Nc-1)

                # Extract the polynomial coefficients
                coefficients[m,:] = V.coeffs

                # Export the data of one reppresentative trajectory
                m==1 && (writeout(hcat(bins, pdf), "../data/distribution/fixed_monostable_saddle_node_$n.csv"))
                m==1 && (writeout(hcat(xs, Vs), "../data/potential/fixed_monostable_saddle_node_$n.csv"))
                m==1 && (writeout([var(u), mean(u), a, b], "../data/settings/$n.csv"))
        end

        # Export the coefficients data
        writeout(coefficients, "../data/coefficients/fixed_monostable_saddle_node_$n.csv")
end

# Execute the additional scripts for this simulation
include("../postprocessing/01_postprocessing.jl")
include("../plotting/01_plotting.jl")
