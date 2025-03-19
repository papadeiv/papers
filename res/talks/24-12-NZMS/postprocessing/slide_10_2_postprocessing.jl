include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

# Read in the timeseries of lattice simulation 
μ = readin("../data/slide_10_2/μ.csv")
ex_traj = readin("../data/slide_10_2/ensemble/1.csv")

# Get the parameters of the ensemble problem
Nμ = length(μ)
Ne = size(ex_traj)[1]

# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e2)
# Number of polynomial coefficients of the least-squares regression
Nc = convert(Int64,4e0)

# Initialise the matrix of polynomial coefficients
coefficients = Matrix{Float64}(undef, Ne, Nc)
mean_c = Matrix{Float64}(undef, Nμ, Nc)

# Loop over the parameter range 
printstyled("Fitting the ensemble data with polynomial least squares\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Read the ensemble solution for the current parameter value
        E = readin("../data/slide_10_2/ensemble/$n.csv")

        # Loop over the particles trajectories
        for m in 1:Ne
                # Extract the trajectory
                local u = E[m,:]

                # Fit an empirical distribution to the data
                bins, pdf = fit_distribution(u, n_bins=Nb+1)

                # Approximate the potential function using the empirical distribution
                xs, Vs = invert_equilibrium_distribution(bins, pdf, 0.125) # The value 0.125 is an artifact and it should be std(u) instead 
                V = approximate_potential(xs, Vs, degree=Nc-1)

                # Extract the polynomial coefficients
                coefficients[m,:] = V.coeffs

                # Export the data of one reppresentative trajectory
                m==1 && (writeout(hcat(xs, Vs), "../data/slide_10_2/potential/$n.csv"))
        end

        # Export the coefficients data
        writeout(coefficients, "../data/slide_10_2/fit/$n.csv")

        # Compute the ensemble mean of the coefficients
        mean_c[n,:] = [mean(coefficients[:,k]) for k in 1:Nc]
end

# Export the ensemble mean of the coefficients
writeout(mean_c, "../data/slide_10_2/fit/mean_coefficients.csv")
