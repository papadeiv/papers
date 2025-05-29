include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/EarlyWarningSignals.jl")

# Define the noise level
σ = 0.200::Float64
D = (σ^2)/(2.0::Float64) 

# Import the data
U_sn = readin("../data/slide_07/saddle_node.csv")
U_OUP = readin("../data/slide_07/OUP.csv")
parameters = readin("../data/slide_07/parameters.csv")
t = parameters[:,1]
μ = parameters[:,2]
ϴ = parameters[:,3]
x0 = parameters[:,4]
Nμ = length(μ)

# Define arrays to store the analytical and sample variance of the 2 solutions
true_var_sn = Vector{Float64}(undef, Nμ)
true_var_OUP = Vector{Float64}(undef, Nμ)
sample_var_sn = Vector{Float64}(undef, Nμ)
sample_var_OUP = Vector{Float64}(undef, Nμ)

# Loop over the parameter values
printstyled("Compute the variance across the sweep using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Nμ
        # Get the timestamp of the last realisation of the solution
        T = convert(Int64,t[n]) - 20
        
        # Extract the solutions at the current parameter value
        u_sn = U_sn[1:T,n] 
        u_OUP = U_OUP[:,n] 

        # Compute the sample variance
        sample_var_sn[n] = var(u_sn)
        sample_var_OUP[n] = var(u_OUP)

        # Compute the analytical variance of the OUP
        true_var_OUP[n] = (σ^2)/(2*abs(ϴ[n]))

        # Define the stationary density of the saddle-node solution
        V(x) =  μ[n]*x + (1.0::Float64/3.0::Float64)*(x^3) - (2.0::Float64/3.0::Float64)*((-μ[n])^(1.5::Float64))
        ρ(x, q) = exp(-(1.0::Float64/D)*(V(x)))
        N = get_normalisation_constant(ρ, (-sqrt(-μ[n]), Inf))
        p(x, q) = N*ρ(x, q)

        # Define the integrand analytical variance of the stationary density above
        f(x, q) = ((x - sqrt(-μ[n]))^2)*p(x, q)

        # Compute the analytical variance of the saddle-node solution
        integral = IntegralProblem(f, (-sqrt(-μ[n]), Inf))
        true_var_sn[n] = (solve(integral, QuadGKJL(); maxiters=10000, reltol=1e-8, abstol=1e-8)).u
end

# Export the data
writeout(hcat(true_var_sn, true_var_OUP, sample_var_sn, sample_var_OUP), "../data/slide_07/variance.csv")
