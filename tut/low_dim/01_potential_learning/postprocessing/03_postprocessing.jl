include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/IO.jl")

# Import the parameters values
μ = readin("../data/parameter.csv")

# Number of fixed values
Nμ = length(μ)-1

# Number of polynomial coefficients of the least-squares fit
d = 4::Int64

# Matrices storing the mean and the variance of the ensemble distribution of the 4 coefficients
coeff_mean = Matrix{Float64}(undef, Nμ, d)
coeff_var = Matrix{Float64}(undef, Nμ, d)

# Loop over the paramater's values
printstyled("Postprocessing\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
          # Read in the coefficients distribution at current parameter value 
        local coefficients = readin("../data/coefficients/fixed_monostable_saddle_node_$n.csv")

        # Loop over the coefficients 
        for m in 1:d
                # Extract the m-th coefficient distribution
                c_m = coefficients[:,m]
                
                # Compute the moments of the distribution
                coeff_mean[n,m] = mean(c_m) 
                coeff_var[n,m] = var(c_m)
        end
end

# Export the ensemble mean and variance of the coefficients
writeout(coeff_mean, "../data/coefficients/mean.csv")
writeout(coeff_var, "../data/coefficients/variance.csv")
