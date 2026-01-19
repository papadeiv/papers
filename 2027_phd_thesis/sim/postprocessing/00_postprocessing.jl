include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the parameter range
μ = readin("../data/parameter.csv")

# Get the number of parameter's values
Nμ = length(μ)

# Loop over the parameter values
printstyled("Postprocessing\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Read all the relevant data and format it
        settings = readin("../data/settings/$n.csv")
        local σ = settings[1]
        local D = (σ^2)/2
        local ϴ = settings[4]

        OUP = readin("../data/distribution/OUP$n.csv")
        bins = OUP[:,1]
        pdf = OUP[:,2]

        histogram = readin("../data/distribution/histogram$n.csv")
        gaussian_bins = histogram[:,1]
        weights = histogram[:,2]

        gaussian_sample = readin("../data/distribution/gaussian$n.csv")
        gaussian_pdf = gaussian_sample[:,2]

        # Get the local potential from the inversion of the stationary distribution formula
        xs, Vs = invert_equilibrium_distribution(bins, pdf, sqrt(D/ϴ))
        hist_xs, hist_Vs = invert_equilibrium_distribution(gaussian_bins, weights, sqrt(D/ϴ))
        gauss_xs, gauss_Vs = invert_equilibrium_distribution(gaussian_bins, gaussian_pdf, sqrt(D/ϴ))

        # Export the data
        writeout(hcat(xs, Vs), "../data/potential/OUP$n.csv")
        writeout(hcat(hist_xs, hist_Vs), "../data/potential/histogram$n.csv")
        writeout(hcat(gauss_xs, gauss_Vs), "../data/potential/gaussian$n.csv")
end
