include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the initial conditions 
x0 = readin("../data/fig:nonautonomous/x0.csv")
Nx = length(x0)

# Import the parameter shift rates 
rates = readin("../data/fig:nonautonomous/rates.csv")
Nε = length(rates)

# Import an example solution
example = readin("../data/fig:nonautonomous/solutions/1.csv")
t = example[:,1]
Nt = length(t)

# Define empty arrays to store postprocessed data 
x3 = Matrix{Float64}(undef, Nt, Nε)
tip_cnt = Vector{Float64}(undef, Nε) 
tip_idx = Matrix{Int64}(undef, Nx, Nε)

# Loop over the rates
printstyled("Postprocessing the data\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε
        # Get the current rate
        ε = rates[n]

        # Import the solutions at the current parameter shift rate 
        solutions = readin("../data/fig:nonautonomous/solutions/$n.csv")
        local t = solutions[:,1]
        local λ = solutions[:,2]
        local X = solutions[:,3:end]

        # Compute the unstable equilibria shift at the current rate
        x3[:,n] = [λ[m] for m in 1:Nt]

        # Initialise count index for the tipped solutions
        cnt = 0::Int64

        # Loop over solutions for the current parameter shift rate
        for m in 1:Nx
                # Extract current solution
                x = X[:,m]

                # Check for irreversible R-tipping
                if x[end] < x3[end,n] 
                        # Update the tipped counter
                        cnt = cnt + 1
                        # Store the index of the tipped solution
                        tip_idx[cnt,n] = m
                end
        end

        # Store the tipped counter at the current parameter shift rate
        tip_cnt[n] = cnt
end

# Export the data
writeout(x3, "../data/fig:nonautonomous/x3.csv")
writeout(tip_cnt, "../data/fig:nonautonomous/tip_cnt.csv")
writeout(tip_idx, "../data/fig:nonautonomous/tip_idx.csv")
