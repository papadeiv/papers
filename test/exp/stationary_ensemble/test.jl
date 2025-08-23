using DifferentialEquations, ForwardDiff, Roots
using LinearAlgebra, Polynomials, Integrals
using StatsBase, LsqFit

# Import functions
include("../../utils/reconstruction_quasipotential.jl")
include("../../utils/analyse_system.jl")
include("../../utils/evolve_system.jl")
include("../../utils/debugging.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")

# Solve the ensemble problem 
t, X = evolve_ensemble(f, η, μ, IC=[x0], δt=δt, Nt=Nt, Ne=Ne)

# Loop over the ensemble's sample paths
printstyled("Computing the least-squares solutions across the ensemble\n"; bold=true, underline=true, color=:light_blue)
#=@showprogress=#for n in 1:Ne
        # Extract the current solution from the ensemble and center it
        u = X[n,:] .- x0 

        # Define a Guassian guess
        guess = [0.0, 1.0, 0.1] + β.*randn(Nc)
        #guess = [ct1, ct2, ct3]
        #guess = β.*randn(Nc)
        debug(n, [""], [guess])

        # Solve the NLLS problem and reconstruct the potential
        coefficients[n,:] = fit_potential(u, n_coeff=Nc, n_bins=Nb, noise=σ, initial_guess=guess)
        c1 = coefficients[n,1] 
        c2 = coefficients[n,2] 
        c3 = coefficients[n,3] 
end
