include("../EscapeProblem.jl")
include("../Plotting.jl")

using CSV, DataFrames
df = DataFrame(CSV.File("../data/analysis.csv"; delim=',', header=false))
μ_range = df[!,1]
escape_rate_analytical = df[!,2]
escape_time_ensemble = df[!,3]
escape_rate_ensemble = Float64[]
all_escaped = Int64[]
last_escape = Float64[]
N_escapes = Float64[]
criteria = Float64[]

df = DataFrame(CSV.File("../data/time.csv"; delim=',', header=false))
t = df[!,1]
T = t[end]

df = DataFrame(CSV.File("../data/ensemble/1.csv"; delim=',', header=false))
Ne = length(df[!,1])+1

treshold = 500::Int64

plot_ensemble(equilibria, μ_range[1], 1)

using ProgressMeter
@showprogress for p in 1:length(μ_range)
        local df = DataFrame(CSV.File("../data/hitting_times/$p.csv"; delim=','))
        h = df[!,1]
        N_exit = length(h)+1
        if N_exit < Ne
                push!(all_escaped, 1::Int64)
        elseif N_exit==Ne
                push!(all_escaped, 0::Int64)
        end
        push!(last_escape, h[end])
        push!(N_escapes, N_exit/Ne)
        escapes = ensemble_escape_rate(T, h, Ne, treshold)
        push!(escape_rate_ensemble, escapes[1])
        push!(criteria, escapes[2])

        # Plot the hitting times distribution
        #plot_escapes(t, h, μ_range[p], p)

        # Get the equilibria at the current parameter value
        f(x, μ) = 2 - 6*x + (12/5)*x^2
        a = [-μ_range[1], -2, +3, -0.8]
        equilibria = get_equilibria(f, a, μ_range[1])

        # Plot the ensemble trajectories
        plot_ensemble(equilibria, μ_range[p], p)
end

#=
using LaTeXStrings, CairoMakie, Makie.Colors
CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)

# Plot on primary y-axis (LEFT)
ax1L = Axis(fig[1,1],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"x",
          xlabelvisible = false,
          ylabelvisible = false,
          xlabelsize = 20,
          ylabelsize = 20,
          xlabelcolor = :black,
          ylabelcolor = :black,
          xlabelpadding = 0.0,
          ylabelpadding = 0.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :left,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
          xticksvisible = true,
          yticksvisible = true,
          xticksize = 6,
          yticksize = 6,
          # Ticks labels
          xticklabelsvisible = true,
          yticklabelsvisible = true,
          xticklabelsize = 18,
          yticklabelsize = 18,
          xticklabelalign = (:right, :top),
          yticklabelalign = (:right, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :blue,
          xtickformat = "{:.6f}",
          ytickformat = "{:.4f}",
)
# Plot the empirical escape rate computed from the ensemble
for n in 1:length(μ_range)
        if criteria[n] == 1::Int64
                scatter!(ax1L, μ_range[n], escape_rate_ensemble[n], color = :lightskyblue1, strokecolor = :black, strokewidth = 1, markersize = 20)
        elseif criteria[n] == 0::Int64
                scatter!(ax1L, μ_range[n], escape_rate_ensemble[n], color = :deepskyblue3, strokecolor = :black, strokewidth = 1, markersize = 20)
        end
end
# Plot the empirical escape rate computed as the inverse of the (mean) escape time of the ensemble
for n in 1:length(μ_range)
        if all_escaped[n] == 1::Int64
                scatter!(ax1L, μ_range[n], 1.0/escape_time_ensemble[n], color = (:lightsalmon,0.95), strokecolor = :black, strokewidth = 1, markersize = 20)
        elseif all_escaped[n] == 0::Int64
                scatter!(ax1L, μ_range[n], 1.0/escape_time_ensemble[n], color = (:firebrick1,0.95), strokecolor = :black, strokewidth = 1, markersize = 20)
        end
end
# Plot the analytical escape rate given by Kramer's formula
lines!(ax1L, μ_range, escape_rate_analytical, color = :blue, linewidth = 1.5)

# Plot on secondary y-axis (RIGHT)
ax1R = Axis(fig[1,1],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), (0.0,750.0)),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"x",
          xlabelvisible = false,
          ylabelvisible = false,
          xlabelsize = 20,
          ylabelsize = 20,
          xlabelcolor = :black,
          ylabelcolor = :black,
          xlabelpadding = 0.0,
          ylabelpadding = 0.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :right,
          # Ticks
          xticks = μ_range,
          #yticks = equilibria[:,1],
          xticksvisible = false,
          yticksvisible = true,
          xticksize = 6,
          yticksize = 6,
          # Ticks labels
          xticklabelsvisible = false,
          yticklabelsvisible = true,
          xticklabelsize = 18,
          yticklabelsize = 18,
          xticklabelalign = (:right, :top),
          yticklabelalign = (:left, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :red,
          xtickformat = "{:.1f}",
          ytickformat = "{:.1f}",
)
hidespines!(ax1R)
hidexdecorations!(ax1R)
# Compute the (mean) escape time from the ensemble
lines!(ax1R, μ_range, escape_time_ensemble, color = :red, linewidth = 1.5)
# Compute the escape time as the inverse of Kramer's escape rate
lines!(ax1R, μ_range, (N_escapes)./(escape_rate_analytical), color = :red, linewidth = 1.5, linestyle = :dash)

# Export the figure
save("../data/$treshold.png", fig)
=#
