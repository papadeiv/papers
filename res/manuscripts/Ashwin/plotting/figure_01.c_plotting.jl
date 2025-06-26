include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

 # Import the data from csv 
solution = readin("../data/figure_01/solutions/1.csv")
t = solution[:,1]
λ = solution[:,2]

ε = readin("../data/figure_01/rates.csv")
Nε = length(ε)

critical_rates = readin("../data/figure_01/critical_rates.csv")
Ic = critical_rates[:,1]
εc = critical_rates[:,2]
Nεc = length(εc)

stable_eq_1 = readin("../data/figure_01/stable_eq_1.csv")
stable_eq_2 = readin("../data/figure_01/stable_eq_2.csv")
unstable_eq = readin("../data/figure_01/unstable_eq.csv")

# Create and customise the state axis 
fig, ax1 = mkfig(size = [1200,900],
                 pad = (60,40,10,30), # Order is: left, right, bottom, top 
                 bg_out = :white,
                 box_position = [1,1],
                 limits = ((λ[1], λ[end]), (-2.1, 2.1)),
                 lab = [L"\mathbf{\lambda}",L"\mathbf{x\,(\lambda(t))}"],
                 lab_size = [30,30],
                 lab_pad = [-30.0,-30.0],
                 x_ticks = [λ[1],λ[end]],
                 y_ticks = [-2,2],
                 ticks_lab_size = [30,30],
                 ticks_lab_xpos = [:right,:top],
                 ticks_lab_trunc = [0,0],
                )

# Create and customise the parameter axis 
fig, ax2 = mkfig(fig = fig,
                 box_position = [2,1],
                 limits = ((t[1], t[end]), (-1.10, 1.10)),
                 lab = [L"\mathbf{t}",L"\mathbf{\lambda(t) = } \textbf{tanh} \mathbf{(\varepsilon\,t)}"],
                 lab_size = [30,30],
                 lab_pad = [-30.0,-30.0],
                 x_ticks = [t[1],t[end]],
                 y_ticks = [-1,1],
                 ticks_lab_size = [30,30],
                 ticks_lab_trunc = [0,0],
                )

# Loop over all the rates
printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε
        # Import the data from csv 
        local solution = readin("../data/figure_01/solutions/$n.csv")
        local t = solution[:,1]
        local μ = solution[:,2]
        local u = solution[:,3]

        # Plot the state's solution 
        lines!(ax1, μ, u, linewidth = 3, color = n, colormap = Makie.Categorical(:Spectral), colorrange = (1,Nε))
        # Plot the parameter's shift
        lines!(ax2, t, μ, linewidth = 4, color = n, colormap = Makie.Categorical(:Spectral), colorrange = (1,Nε))
end

# Plot the bifurcation diagram of the frozen system
lines!(ax1, stable_eq_1[:,1], stable_eq_1[:,2], linewidth = 4, color = :black)
lines!(ax1, stable_eq_2[:,1], stable_eq_2[:,2], linewidth = 4, color = :black)
lines!(ax1, unstable_eq[:,1], unstable_eq[:,2], linewidth = 4, color = :black, linestyle = :dash)

# Plot the colorbar for all the rates
Colorbar(fig[1:2,2],
         size = 45, 
         spinewidth = 5.0,
         limits = (ε[1], ε[end]), 
         colormap = cgrad(:Spectral, Nε, categorical = true),
         ticks = ε,
         ticksize = 22,
         tickwidth = 5.0,
         tickformat = "{:.$(2)f}",
         ticklabelsize = 30,
         label = L"\mathbf{\varepsilon}",
         labelpadding = 0.0,
         labelsize = 30
        )

# Export the figure 
save("../fig/figure.png", fig)
