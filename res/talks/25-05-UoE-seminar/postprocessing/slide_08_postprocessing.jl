include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/EarlyWarningSignals.jl")

# Define the noise level
σ = 0.250::Float64
D = (σ^2)/(2.0::Float64) 

# Import the data
solution = readin("../data/slide_08/solution.csv")
μ = solution[:,2]
u = solution[:,3]
Nμ = length(μ)

# Define array to store the normalisation constant of the stationary distribution 
N = Vector{Float64}(undef, Nμ)
# Define array to store the probability to escape
tip_prob = Vector{Float64}(undef, Nμ)

# Loop over the parameter values
printstyled("Compute the normalisation constant across using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Nμ
        # Define the analytic potential function
        V(x) = μ[n]*x + x^2 - x^3 + (1/5)*(x^4)

        # Define the unnormalised equilibrium density
        ρ(x, q) = exp(-(1.0::Float64/D)*(V(x)))

        # Compute the normalisation constant
        N[n] = get_normalisation_constant(ρ, (-Inf, Inf))

        # Define the potential's second derivative
        Vxx(x) = 2 - 6*x + (12/5)*(x^2) 

        # Define the tipping boundaries
        f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
        eq = get_equilibria(f, μ[n], domain=[-20,20])
        b = sort(eq[1], rev=true)[1]
        a = 0.0::Float64
        if length(eq[2]) > 0
                a = eq[2][1]
        else
                a = 2.12::Float64 
        end

        # Compute the tipping probability
        if b > a
                estimate = kramer_escape(V, Vxx, a, b, σ)
                if estimate > 1.0
                        tip_prob[n] = 1.0
                else
                        tip_prob[n] = kramer_escape(V, Vxx, a, b, σ)
                end
        else
                tip_prob[n] = 1.0 
        end
end

# Export the data
writeout(hcat(N, tip_prob), "../data/slide_08/escapes.csv")
