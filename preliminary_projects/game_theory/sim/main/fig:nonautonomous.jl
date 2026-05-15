include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#=
# Define the dynamical system (drift of the state variable) 
f(x, λ) = x*(1.0 - x)*(x - λ)

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define a set of initial conditions (ICs) 
Nx = 35
ICs = collect(LinRange(0.05, 0.95, Nx))

# Define the the past limit system
T = 10.0
t∞ = [-T, +T]
Nt = convert(Int64,1e4)

# Define the parameter's rates
Nε = 1000
rates = collect(LinRange(0.01, 0.30, Nε)) 

# Loop over the rates
printstyled("Solving the non-autonomous ODEs for different rates\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε 
        # Define the dynamical system (drift of the parameter)
        k = 0.5
        δ = 0.0
        ε = rates[n]
        g(t) = k*ε*(sech(ε*t + δ))^2

        # Define empty matrix to store the solutions
        X = Matrix{Float64}(undef, (Nt+1), (Nx+2))

        # Loop over different ICs
        for m in 1:length(ICs)
                # Extract the current initial condition
                λ∞ = 0.00247
                x∞ = ICs[m]
                u∞ = [x∞, λ∞]

                # Solve the non-autonomous ODE for the current initial condition
                local t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=Nt)

                # Store the results in the matrix
                if m==1
                        X[:,1] = t
                        X[:,2] = λ
                end
                X[:,m+2] = x
        end

        # Export the data 
        writeout(X, "../data/fig:nonautonomous/solutions/$n.csv")
end

# Export the initial conditions and parameter rates
writeout(ICs, "../data/fig:nonautonomous/x0.csv")
writeout(rates, "../data/fig:nonautonomous/rates.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/fig:nonautonomous.jl")
=#
include("../plotting/fig:nonautonomous.jl")
