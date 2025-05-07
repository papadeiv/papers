include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Import the timestamps and parameter values
time = readin("../data/figure_05/time.csv")
Nt = length(time)
solutions = readin("../data/figure_05/solutions.csv")
Ne = length(solutions[:,1])

# Define the time index of the truncation of the timeseries
Ns = 850

# Vector to store the escape time
T_esc = Float64[]
E_esc = Int64[]
# Vector to store the tipping time
T_tip = Float64[]
E_tip = Int64[]
# Vector to store the trajectories whose escape time is close to the ensemble mean
E_mean = Int64[]

# Loop over the ensemble
for e in 1:Ne
        # Compute the time of escape
        check, escape = check_escapes(solutions[e,1:Ns], [-100,0.867])

        # Check if the trajectory has escapes
        if check
                # Update the escape time vectors
                push!(T_esc, time[escape])
                push!(E_esc, e)
        end

        # Compute the time it takes to settle onto the alternative equilibrium
        check, tipping = check_escapes(solutions[e,1:Ns], [-100,2.8])

        # Check if the trajectory has tipped 
        if check
                # Check if the trajectory has not returned to the basin
                returned, escape = check_escapes(solutions[e,Ns:end], [0.867,+100])
                if !returned
                        # Update the escape time vectors
                        push!(T_tip, time[tipping])
                        push!(E_tip, e)
                end
        end
end

# Fit an histogram to the escape time distribution
bins, pdf = fit_distribution(T_esc, n_bins=41)

# Compute the ensemble mean escape time
T_mean = mean(T_esc)

# Extract the trajectories whose escape time is within some distance from the mean
for e in 1:length(E_esc)
        if T_esc[e] > (0.99*T_mean) && T_esc[e] < (1.01*T_mean)
                check = findall(i -> i==E_esc[e], E_tip)
                if length(check) > 0
                        push!(E_mean, E_esc[e])
                end
        end
end

# Export the data
writeout(hcat(E_esc, T_esc), "../data/figure_05/escapes.csv")
writeout(hcat(E_tip, T_tip), "../data/figure_05/tipping.csv")
writeout(E_mean, "../data/figure_05/mean_escaped.csv")
writeout(hcat(bins, pdf), "../data/figure_05/escape_time_distribution.csv")
