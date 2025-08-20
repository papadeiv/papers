using CairoMakie, Makie.Colors
using CSV, DataFrames

CairoMakie.activate!(; px_per_unit = 2)

df = DataFrame(CSV.File("./Yang23.csv"; delim=','))
time = df[!,1]
SST = df[!,2]

fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
ax = Axis(fig[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    #limits = ((0,10),(18,27)),
    # Title
    title = "heavy Fermi phase transition",
    titlevisible = true,
    titlesize = 50,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = "temperature (K)",
    xlabelvisible = true,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -60.0,
    xticks = [2,200],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 50,
    xtickformat = "{:.0f}",
    xscale = log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = "resonance (THz)",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    yticks = [0.25,3.75],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.2f}",
    yscale = identity,
    yaxisposition = :left,
)
#text!(3, 0.5, text = L"\text{(e)}", align = [:right, :top], color = :black, fontsize = 70)
#text!(9, 1, text = L"\times 10 ^{-3}", align = [:right, :bottom], color = :black, fontsize = 30)
lines!(ax, time, SST, color = :snow4, linewidth = 3)
save("../results/physics/Yang23.png", fig)

