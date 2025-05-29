include("../../../../inc/IO.jl")
include("../../../../inc/LatticeDynamics.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Import the timestamps and parameter values
time = readin("../data/figure_05/time.csv")
Nt = length(time)
parameter = readin("../data/figure_05/parameter.csv")
Nμ = length(parameter)

#=
# Array to store the various metrics of the solution
L2norm = Vector{Float64}(undef, Nμ)
u_mean = Vector{Float64}(undef, Nμ)
u = Vector{Float64}(undef, Nμ)

# Loop over the parameter sweep
printstyled("Postprocessing the data\n"; bold=true, underline=true, color=:green)
@showprogress for n in 1:Nμ
        # Import the solution
        solution = readin("../data/figure_05/solutions/$n.csv") 

        # Compute the L2 norm of the steady-state
        L2norm[n] = norm(solution[:,end], 2)  

        # Compute the spatial mean of the steady-state
        u_mean[n] = mean(solution[:,end])

        # Update the value at the midpoint of the spatial domain at steady-state
        u[n] = solution[500,end]
end
=#

metrics = readin("../data/figure_05/metrics.csv")
# Compute the variance of the spatial mean
t_var, u_var = ews.variance(parameter, metrics[:,1], 0.05)

# Export the data
#writeout(hcat(L2norm, u_mean, u), "../data/figure_05/metrics.csv")
writeout(hcat(t_var, u_var), "../data/figure_05/variance.csv")
