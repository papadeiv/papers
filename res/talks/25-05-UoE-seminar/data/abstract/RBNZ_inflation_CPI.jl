using CairoMakie, Makie.Colors
using CSV, DataFrames

df = DataFrame(CSV.File("./RBNZ_inflation_CPI.csv"; delim=',', types=Float64, header=false))
CPI = df[!,3]

CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)
ax = Axis(fig[1,1],
           # Background
           #backgroundcolor = :transparent,
           xgridvisible = false,
           ygridvisible = false,
           limits = (nothing, nothing),
           # Title
           #title = L"\mu = %$parameter",
           titlevisible = false,
           titlesize = 22,
           titlealign = :center,
           titlegap = 4.0,
           # Axes labels
           xlabel = "time",
           ylabel = "HPI",
           xlabelvisible = true,
           ylabelvisible = true,
           xlabelsize = 20,
           ylabelsize = 20,
           xlabelcolor = :black,
           ylabelcolor = :black,
           xlabelpadding = 0.0,
           ylabelpadding = 0.0,
           # Axes scale, position and direction
           xscale = identity, 
           yscale = identity, #log10,
           xreversed = false,
           yreversed = false,
           xaxisposition = :bottom,
           yaxisposition = :left,
           # Ticks
           #xticks = LinRange(μ0,-0.371231,10),
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
           yticklabelcolor = :black,
           xtickformat = "{:.2f}",
           ytickformat = "{:.4f}",
)
lines!(ax, LinRange(1,length(CPI),length(CPI)), CPI, linewidth = 3)
save("./fig.png", fig)

