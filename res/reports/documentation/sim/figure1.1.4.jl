using CairoMakie, Makie.Colors
using CSV, DataFrames

CairoMakie.activate!(; px_per_unit = 2)

df = DataFrame(CSV.File("./Dai12.csv"; delim=','))
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
    title = "lab-specie collapse",
    titlevisible = true,
    titlesize = 50,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = "days",
    xlabelvisible = true,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -60.0,
    xticks = [2,3,4,7,8,9],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 50,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = "population density",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -60.0,
    yticks = [0.04,1],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    ytickformat = "{:.2f}",
    yscale = log10,
    yaxisposition = :left,
)
#text!(2.5, 0.06, text = L"\text{(d)}", align = [:right, :top], color = :black, fontsize = 70)
#text!(9, 1, text = L"\times 10 ^{-3}", align = [:right, :bottom], color = :black, fontsize = 30)
lines!(ax, time, SST, color = :orange, linewidth = 3)
save("../results/ecosystems/Dai12.png", fig)

