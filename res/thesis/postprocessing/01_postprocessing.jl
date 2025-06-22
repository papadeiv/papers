include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/IO.jl")

# Import the parameters values
μ = readin("../data/parameter.csv")

# Number of fixed values
Nμ = length(μ)
# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e2)
# Number of polynomial coefficients of the least-squares regression
Nc = convert(Int64,4e0)

# Loop over the paramater's values
printstyled("Postprocessing\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        ####################### 
        # Stationary solution #
        ####################### 
        
        # Import the data
        solution = readin("../data/solutions/monostable_$n.csv")
        # Fit an empirical distribution
        bins, pdf = fit_distribution(solution[:,2], n_bins=Nb+1)
        # Invert the equilibrium distribution
        xs, Vs = invert_equilibrium_distribution(bins, pdf, std(solution[:,2]))
        # Derive a polynomial approximation through linear least-squares
        V = approximate_potential(xs, Vs, degree=Nc-1)
        # Export the data
        writeout(hcat(bins, pdf), "../data/distribution/monostable_$n.csv")
        writeout(hcat(xs, Vs), "../data/potential/monostable_$n.csv")
        writeout(V.coeffs, "../data/coefficients/monostable_$n.csv")

        ######################
        # Detrended solution #
        ######################
        
        # Detrend the solution
        detrend(solution[:,2], "../data/solutions/detrended_$n.csv") 
        residuals = readin("../data/solutions/detrended_$n.csv")
        # Fit an empirical distribution 
        bins, pdf = fit_distribution(residuals, n_bins=Nb+1)
        # Export the data
        writeout(hcat(bins, pdf), "../data/distribution/detrended_monostable_$n.csv")

        ################################
        # Topologically equivalent OUP #
        ################################

        # Import the data
        OUP = readin("../data/solutions/OUP_$n.csv")
        # Fit an empirical distribution
        bins, pdf = fit_distribution(OUP[:,2], n_bins=Nb+1)
        # Export the data
        writeout(hcat(bins, pdf), "../data/distribution/equivalent_OUP_$n.csv")
end
