include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Payoff matrix
A = [[1.0, 0.5, 0.0, 1.0], # Coordination game
     [0.0, 2.0, 1.0, 3.0], # Dominant strategy (x1 stable)
     [1.0, 3.0, 0.0, 2.0], # Dominant strategy (x2 stable)
     [0.5, 1.0, 1.0, 0.0]] # Anti-coordination
Nμ = length(A)

# Define a set of initial conditions (ICs) 
Nx = 100 
ICs = collect(LinRange(0.05, 0.95, Nx))

# Define time parameters of the simulation
Nt = 1000
dt = 1e-2
T = Nt*dt

# Loop over the parameter values
printstyled("Solving the ODE for different payoffs\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Extract the current parameter vector
        μ = A[n]

        # Define the dynamical system (deterministic drift)
        f(x, μ, t) = x*(1 - x)*((μ[1] + μ[4] - μ[2] - μ[3])*x + μ[2] - μ[4])

        # Define empty matrix to store the solutions
        X = Matrix{Float64}(undef, (Nt+1), (Nx+1))

        # Loop over different ICs
        for m in 1:length(ICs)
                # Extract the current initial condition
                x0 = ICs[m]
                
                # Solve the differential equation for the current IC
                dynamics = ODEProblem(f, x0, (0.0, T), μ)
                sol = solve(dynamics, dt=dt, saveat=dt)

                # Extract the solution
                t = sol.t
                x = [(sol.u)[t][1] for t in 1:length(t)]

                # Store the results in the matrix
                if m==1
                        X[:,1] = t
                end
                X[:,m+1] = x
        end

        # Export the solution matrix 
        writeout(X, "../data/figure_01/$n.csv")
end

# Export the ICs
writeout(ICs, "../data/figure_01/x0.csv")

# Execute the plotting script
include("../plotting/figure_01_plotting.jl")
