using CairoMakie, Makie.Colors
using CSV, DataFrames
using Statistics

df = DataFrame(CSV.File("./RBNZ_housing_HPI.csv"; delim=',', groupmark=',', types=Int64, header=false))
consents = df[!,1]
HPI = df[!,2]

variance = Float64[]
window = 10
for n in 1:(length(HPI)-window)
        push!(variance, var(HPI[n:(n+window)]))
end

CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
axL = Axis(fig[1,1],
           # Background
           #backgroundcolor = :transparent,
           xgridvisible = false,
           ygridvisible = false,
           limits = ((1,138), nothing),
           # Title
           #title = L"\mu = %$parameter",
           titlevisible = false,
           titlesize = 22,
           titlealign = :center,
           titlegap = 4.0,
           # Axes labels
           xlabel = "time",
           ylabel = "House price index (HPI)",
           xlabelvisible = true,
           ylabelvisible = true,
           xlabelsize = 20,
           ylabelsize = 20,
           xlabelcolor = :black,
           ylabelcolor = :dodgerblue4,
           xlabelpadding = 0.0,
           ylabelpadding = 15.0,
           # Axes scale, position and direction
           xscale = identity, 
           yscale = identity, #log10,
           xreversed = false,
           yreversed = false,
           xaxisposition = :bottom,
           yaxisposition = :left,
           # Ticks
           xticks = ([8,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136],
                     ["1991","1993","1995","1997","1999","2001","2003","2005","2007","2009","2011","2013","2015","2017","2019","2021","2023"]),
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
           xticklabelrotation = pi/4,
           yticklabelrotation = 0,
           xticklabelcolor = :black,
           yticklabelcolor = :dodgerblue4,
           #xtickformat = "{:.0f}",
           ytickformat = "{:.0f}",
)
lines!(axL, LinRange(1,length(HPI),length(HPI)), HPI, color = :dodgerblue4, linewidth = 3)
#=
axR = Axis(fig[1,1],
           # Background
           #backgroundcolor = :transparent,
           xgridvisible = false,
           ygridvisible = false,
           limits = ((1,138), nothing),
           # Title
           #title = L"\mu = %$parameter",
           titlevisible = false,
           titlesize = 22,
           titlealign = :center,
           titlegap = 4.0,
           # Axes labels
           xlabel = "time",
           ylabel = "Building consents",
           xlabelvisible = true,
           ylabelvisible = true,
           xlabelsize = 20,
           ylabelsize = 20,
           xlabelcolor = :black,
           ylabelcolor = :brown2,
           xlabelpadding = 0.0,
           ylabelpadding = 0.0,
           # Axes scale, position and direction
           xscale = identity, 
           yscale = identity, #log10,
           xreversed = false,
           yreversed = false,
           xaxisposition = :bottom,
           yaxisposition = :right,
           # Ticks
           #xticks = LinRange(μ0,-0.371231,10),
           #yticks = equilibria[:,1],
           xticksvisible = true,
           yticksvisible = true,
           xtickcolor = :transparent,
           ytickcolor = :black,
           xticksize = 6,
           yticksize = 6,
           # Ticks labels
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xticklabelsize = 18,
           yticklabelsize = 18,
           xticklabelalign = (:right, :top),
           yticklabelalign = (:left, :center),
           xticklabelcolor = :black,
           yticklabelcolor = :brown2,
           xtickformat = "{:.0f}",
           ytickformat = "{:.0f}",
)
hidespines!(axR)
hidexdecorations!(axR)
lines!(axR, 21:1:138, consents[21:end], color = (:brown2,0.99), linewidth = 3)
=#
axR = Axis(fig[1,1],
           # Background
           #backgroundcolor = :transparent,
           xgridvisible = false,
           ygridvisible = false,
           limits = ((1,138), nothing),
           # Title
           #title = L"\mu = %$parameter",
           titlevisible = false,
           titlesize = 22,
           titlealign = :center,
           titlegap = 4.0,
           # Axes labels
           xlabel = "time",
           ylabel = "Variance",
           xlabelvisible = true,
           ylabelvisible = true,
           xlabelsize = 20,
           ylabelsize = 20,
           xlabelcolor = :black,
           ylabelcolor = :brown2,
           xlabelpadding = 0.0,
           ylabelpadding = 0.0,
           # Axes scale, position and direction
           xscale = identity, 
           yscale = identity, #log10,
           xreversed = false,
           yreversed = false,
           xaxisposition = :bottom,
           yaxisposition = :right,
           # Ticks
           #xticks = LinRange(μ0,-0.371231,10),
           #yticks = equilibria[:,1],
           xticksvisible = true,
           yticksvisible = true,
           xtickcolor = :transparent,
           ytickcolor = :black,
           xticksize = 6,
           yticksize = 6,
           # Ticks labels
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xticklabelsize = 18,
           yticklabelsize = 18,
           xticklabelalign = (:right, :top),
           yticklabelalign = (:left, :center),
           xticklabelcolor = :black,
           yticklabelcolor = :brown2,
           xtickformat = "{:.0f}",
           ytickformat = "{:.0f}",
)
hidespines!(axR)
hidexdecorations!(axR)
lines!(axR, (window+1):1:138, variance, color = (:brown2,1.00), linewidth = 3)
save("./fig.png", fig)
