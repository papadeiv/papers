include("../EscapeProblem.jl")

σ=0.010

μ0 = -0.3650 #-0.371
μ_range=LinRange(μ0,-0.371231,200)

using CSV, Tables
mat = Matrix{Float64}(undef, length(μ_range), 3)
mat[:,1] = μ_range

escape_rate_analytical = Float64[] 
escape_rate_ensemble = Float64[] 
escape_time_analytical = Float64[] 
escape_time_ensemble = Float64[] 
last_escape_time = Int64[]

using ProgressMeter, Statistics
@showprogress for p in 1:length(μ_range)
        μ = μ_range[p]

        equilibria = get_equilibria(μ)

        push!(escape_rate_analytical, Kramer(equilibria[1,1], equilibria[2,1], μ, σ))

        t, u = evolve_ensemble(μ, σ, equilibria)
        if p==1
                CSV.write("./data/time.csv", Tables.table(t), delim=',', writeheader=false)
        end
        CSV.write("./data/ensemble/$p.csv", Tables.table(u), delim=',', writeheader=false)

        hitting_times = ensemble_escapes(t, u, equilibria)
        CSV.write("./data/hitting_times/$p.csv", Tables.table(hitting_times), delim=',', writeheader=false)

        if length(hitting_times) > 0
                push!(escape_time_ensemble, mean(hitting_times))
        else
                push!(escape_time_ensemble, t[end])
        end
end

mat[:,2] = escape_rate_analytical 
mat[:,3] = escape_time_ensemble
CSV.write("./data/analysis.csv", Tables.table(mat), delim=',', writeheader=false)
