using CairoMakie, Makie.Colors
using CSV, DataFrames

CairoMakie.activate!(; px_per_unit = 2)

df = DataFrame(CSV.File("./Tripati05.csv"; delim=','))
time = df[!,1]
SST = df[!,2]

fig = Figure(; size = (1200, 1200), backgroundcolor = :transparent)
ax = Axis(fig[1,1:4],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    #limits = ((0,10),(18,27)),
    # Title
    title = "end of Greenhouse Earth",
    titlevisible = true,
    titlesize = 50,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = "Myears BCE",
    xlabelvisible = true,
    xlabelsize = 50,
    xlabelcolor = :black,
    xlabelpadding = -60.0,
    xticks = [39,37,34,32],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 50,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = "Carbonate",
    ylabelvisible = true,
    ylabelsize = 50,
    ylabelcolor = :black,
    ylabelpadding = -48.0,
    yticks = [0,15,70,85],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 50,
    #ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
ax.xreversed = true
#text!(32.5, 14, text = L"\text{(a)}", align = [:left, :top], color = :black, fontsize = 70)
#text!(9, 1, text = L"\times 10 ^{-3}", align = [:right, :bottom], color = :black, fontsize = 30)
lines!(ax, time, SST, color = :blue, linewidth = 3)
save("../../results/climate/paleoclimate/Tripati05.png", fig)

