include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the timeseries data 
variance = readin("../data/slide_07/variance.csv")
parameters = readin("../data/slide_07/parameters.csv")
t = parameters[:,1]
μ = parameters[:,2]
ϴ = parameters[:,3]
x0 = parameters[:,4]
Nμ = length(μ)

printstyled("Generating figures using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
# Create and customise the OUP variance axis 
fig, ax = mkfig(size = [1600,800],
                bg_out = colorant"rgb(147,237,213)",
                box_position = [1,1],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((μ[1], 0), (0, 0.1)),
                title = L"\textbf{OUP}",
                toggle_title = true,
                title_color = RGBAf(0,0.236,0.236),
                lab = [L"\mathbf{\mu}", L"\textbf{Var}"],
                lab_pad = [-60.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                x_ticks = [μ[1],0],
                y_ticks = [0,0.1],
                ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                ticks_lab_trunc = [0,1]
               )
# Plot the OUP sample variance 
lines!(ax, μ, variance[:,4], linewidth = 1, color = (colorant"rgb(38,106,176)", 1.0))
# Plot the OUP analytic variance 
lines!(ax, μ, variance[:,2], linewidth = 4.5, color = (:teal, 1.0))

# Create and customise the saddle-node variance axis 
fig, ax = mkfig(fig = fig,
                #bg_out = colorant"rgb(147,237,213)",
                box_position = [1,2],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((μ[1], 0), (0, 0.1)),
                title = L"\textbf{saddle-node}",
                toggle_title = true,
                title_color = RGBAf(0,0.236,0.236),
                lab = [L"\mathbf{\mu}", L"\textbf{Var}"],
                toggle_lab = [true,false],
                lab_pad = [-60.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                x_ticks = [μ[1],0],
                y_ticks = [0,0.1],
                ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                toggle_ticks_lab = [true,false],
                ticks_lab_trunc = [0,1]
               )
# Plot the saddle-node sample variance 
lines!(ax, μ, variance[:,3], linewidth = 1, color = (colorant"rgb(38,106,176)", 1.0))
# Plot the saddle-node analytic variance 
lines!(ax, μ, variance[:,1], linewidth = 4.5, color = (:teal, 1.0))

# Export the figure
save("../fig/slide_07.png", fig)
