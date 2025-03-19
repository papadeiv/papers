include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")

# Read in the data of the paramater values 
μ = readin("../data/slide_10_1/μ.csv")

# Get the parameters of the ensemble problem
Nμ = length(μ)

# Number of bins for the histogram of the realizations
Nb = convert(Int64,1e3)

# Loop over the parameter range 
printstyled("Fitting an empirical distribution to the solutions of the SDEs\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Read the solution for the current parameter value
        local u = readin("../data/slide_10_1/solutions/$n.csv")

        # Fit an empirical distribution to the data
        bins, pdf = fit_distribution(u, interval=[2,4], n_bins=Nb+1)

        # Export the data
        writeout(hcat(bins, pdf), "../data/slide_10_1/distribution/$n.csv")
end
