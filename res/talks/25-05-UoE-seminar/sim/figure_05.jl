include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/LatticeDynamics.jl")

# Spatial domain
L = 2*π 
Nx = convert(Int64, 1e3)
dx = L/Nx

# Temporal domain
T = 100.0::Float64
δt = 1e-1 

# Parameter sweep
μ0 = -(0.50::Float64)
μf = 0.50::Float64
Nμ = convert(Int64,1e3)
μ = LinRange(μ0, μf, Nμ)

# Discretised spatial differential operators
∂2x = Tridiagonal([1.0::Float64 for n in 1:(Nx-1)], [-(2.0::Float64) for n in 1:Nx], [1.0::Float64 for n in 1:(Nx-1)]) 
∂4x = BandedMatrix{Float64}(undef, (Nx,Nx), (2,2))
∂4x[band(0)] .= 6.0::Float64
∂4x[band(1)] .= ∂4x[band(-1)] .= -(4.0::Float64)
∂4x[band(2)] .= ∂4x[band(-2)] .= 1.0::Float64

# Initial condition
σ = 0.050::Float64
u0 = 0.5::Float64.*ones(Float64, Nx) .+ σ.*rand(Float64, Nx)

# Boundary conditions
∂2x[2,1] = ∂2x[end-1,end] = 2.0::Float64
∂4x[2,1] = ∂4x[end-1,end] = -(6.0::Float64) 
∂4x[3,1] = ∂4x[end-2,end] = 4.0::Float64 

# Drift
ɑ = 2.0::Float64
β = 1.0::Float64
q = 1.0::Float64
f(u, μ, j) = -2*dot(∂2x[j,:],u) - dot(∂4x[j,:],u) + (μ - q^2)*u[j] +ɑ*(u[j])^3 - β*(u[j])^5

# Diffusion
η(u, μ, j) = σ

# Perform the parameter sweep
printstyled("Performing the parameter sweep\n"; bold=true, underline=true, color=:green)
@showprogress for n in 1:Nμ
        # Define the functions
        function drift!(dudt, u, p, t)
                for j in 1:Nx
                        dudt[j] = f(u, μ[n], j)
                end
        end
        function diffusion!(dudt, u, p, t)
                for j in 1:Nx
                        dudt[j] = η(u, μ[n], j)
                end
        end

        # Define the SDE 
        sde = SDEProblem(drift!, diffusion!, u0, (0.0, T))

        # Propagate the SDE forward in time
        sol = solve(sde, EM(), dt=δt, verbose=false)#, saveat=1000*δt)

        # Extract the timestamps
        time = sol.t
        local Nt = length(time)

        # Extract the solution
        state = Matrix{Float64}(undef, Nx, Nt)
        for t in 1:Nt
                state[:,t] = (sol.u)[t]
        end

        # Export the data
        n==1 && (writeout(time, "../data/figure_05/time.csv")) 
        n==1 && (writeout(μ, "../data/figure_05/parameter.csv")) 
        writeout(state, "../data/figure_05/solutions/$n.csv")
end

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_05_postprocessing.jl")
include("../plotting/figure_05_plotting.jl")
