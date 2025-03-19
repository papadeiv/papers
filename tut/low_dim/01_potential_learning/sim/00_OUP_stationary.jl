include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/IO.jl")

# Number of realizations 
Nt = convert(Int64,1e5)
# Number of (fixed) parameter values
Nμ = convert(Int64,1e1)
# Number of bins for the histogram of the realizations
Nb = convert(Int64,2e1)
# Polynomial degree of the least-squares regression 
Nc = convert(Int64,2e0)

# Define the parameter of the OUP
ϴ = 0.1::Float64
# Define the dynamics of an OUP
f(x, μ) = -ϴ*(x - μ)

# Specify the (additive) noise level in the OUP
σ = 0.125::Float64 
# Rdefine the noise in terms of the (stochastic) diffusion coefficient 
D = σ^2/2
g(x) = sqrt(2*D) 

# Specify the parameter range and export it 
μ = LinRange(0.000,1.000,Nμ)
writeout(μ, "../data/parameter.csv")

# Initialise the matrix of polynomial coefficients
coefficients = Matrix{Float64}(undef, Nμ, Nc)

# Loop over the parameter range 
printstyled("Simulating the SDE (OUP) and generating gaussian samples\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Propagate the OUP forward in time
        t, u = evolve_forward_1d(Ux, g, μ[n], Nt=Nt)

        # Fit an empirical distribution to the data
        bins, pdf = fit_distribution(u, n_bins=Nb+1)
        
        # Export the data
        writeout(hcat(bins, pdf), "../data/distribution/OUP$n.csv")

        # Here you want to show that the distribution of a stationary OUP process is equivalent to a gaussian distribution with prescribed mean and variance 
        α = sqrt(D/ϴ) # Variance of the gaussian samples 
        β = μ[n] # Mean of the gaussian samples
        samples = α.*randn(Nt) .+ β # Transformation of a normally distributed random variable with mean 0 and variance 1

        # Extrema of the samples
        a = minimum(samples)
        b = maximum(samples)

        # Range of the samples
        I = b - a

        # Define the left edges of the bins of the histogram
        x = LinRange(a - 0.05*I, b + 0.05*I, Nb+1)

        # Compute the length of each bin (they're all equally spaced)
        dx = abs(x[2]-x[1]).*ones(Nb)

        # Compute the midpoints of the bins (for plotting purposes) 
        gaussian_bins = [(x[n+1]+x[n])/2 for n in 1:Nb]

        # Fit an histogram
        hist = StatsBase.fit(Histogram, samples, x)

        # The weights of the histogram represent the counts of the samples in each bin
        weights = hist.weights

        # Compute the integral of the histogram (normalisation constant)
        N = LinearAlgebra.dot(dx, weights)

        # Safety check
        if N ≈ LinearAlgebra.norm(hist)
                # Export the data
                writeout(hcat(gaussian_bins, weights), "../data/distribution/histogram$n.csv")
        else
                break
        end

        # Take the weights of the histogram and divide them by the norm of said histogram
        gaussian_pdf = weights./N

        # Compute the integral of the normalised distribution above (should be 1!)
        L1_norm = LinearAlgebra.dot(dx, gaussian_pdf)

        # Safety check
        if L1_norm ≈ LinearAlgebra.norm(LinearAlgebra.normalize(hist, mode = :pdf)) && L1_norm ≈ 1.0::Float64
                # Export the data
                writeout(hcat(gaussian_bins, gaussian_pdf), "../data/distribution/gaussian$n.csv")
        else
                break
        end
        
        # Export the settings
        writeout([σ, a, b, ϴ], "../data/settings/$n.csv")
end

# Execute the additional scripts for this simulation
include("../postprocessing/00_postprocessing.jl")
include("../plotting/00_plotting.jl")
