include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

# Payoff matrix
A = [[1.0, 0.5, 0.0, 1.0], # Coordination game
     [0.0, 2.0, 1.0, 3.0], # Dominant strategy (x1 stable)
     [1.0, 3.0, 0.0, 2.0], # Dominant strategy (x2 stable)
     [0.5, 1.0, 1.0, 0.0]] # Anti-coordination
Nμ = length(A)

# Extract x3
x3(μ) = (μ[4] - μ[2])/(μ[1] + μ[4] - μ[2] - μ[3])

# Import the initial conditions 
x0 = readin("../data/fig:autonomous/x0.csv")
Nx = length(x0)

#######################
#   Bifurcation set   #
#######################

# Define limits of the set
ɑ_min = -(2.0::Float64)
ɑ_max = 2.0::Float64
β_min = -(2.0::Float64)
β_max = 2.0::Float64

# Define the ɑ-domains for the admissible set Γ*
ɑ_dom_1 = LinRange(ɑ_min, 0, 1000)
ɑ_dom_2 = LinRange(0, ɑ_max, 1000)

# Define symmetric position for printing the text
δx = 1.0
δy = 1.5

# Create and customise the figure 
fig, ax = mkfig(size = [1000,1000],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((ɑ_min,ɑ_max), (β_min,β_max)),
                lab = [L"\mathbf{\alpha := a - c}",L"\mathbf{\beta := d - b}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [ɑ_min,ɑ_max],
                y_ticks = [β_min,β_max],
                ticks_lab_trunc = [0,0],
               )
# Shade the admissible set Γ*
band!(ɑ_dom_1, ones(1000).*β_min, zeros(1000), color = (:gray, 0.5))
band!(ɑ_dom_2, zeros(1000), ones(1000).*β_max, color = (:gray, 0.5))
# Plot the separatrices of the bifurcation set 
lines!(ax, [0,0], [β_min,β_max], linewidth = 6, color = :black)
lines!(ax, [ɑ_min,ɑ_max], [0,0], linewidth = 6, color = :black)
# Print the name of the regions in the bifurcation set
text!(δx, δy, text = L"\textbf{coordination game}", fontsize = 25, align = (:center,:baseline))
text!(-δx, δy, text = L"\textbf{dominant strategy (}\mathbf{x_1}\textbf{ stable)}", fontsize = 25, align = (:center,:baseline))
text!(δx, -δy, text = L"\textbf{dominant strategy (}\mathbf{x_2}\textbf{ stable)}", fontsize = 25, align = (:center,:baseline))
text!(-δx, -δy, text = L"\textbf{anti-coordination}", fontsize = 25, align = (:center,:baseline))
# Loop over the different strategies
for n in 1:Nμ
        # Compute the position in the bifurcation set
        ɑ = A[n][1] - A[n][3] 
        β = A[n][4] - A[n][2]

        # Plot the position in the bifurcation set
        scatter!(ax, ɑ, β, color = :purple, markersize = 30, strokecolor = :black, strokewidth = 2)
end

# Export the figure 
save("../fig/fig:autonomous_bifset.png", fig)

#########################
#   Coordination game   #
#########################

# Import the solution data at current timestep from csv 
solution = readin("../data/fig:autonomous/1.csv")
t = solution[:,1]
X = solution[:,2:end]

# Create and customise the figure 
fig, ax = mkfig(size = [1000,1000],
                bg_out = :white,
                pad = (30,60,10,5), # Order is: left, right, bottom, top 
                box_position = [1,2],
                limits = ((t[1], t[end]), (-0.05, 1.05)),
                title = L"\textbf{Coordination game}",
                toggle_title = true,
                title_size = 30,
                toggle_lab = [false,false],
                x_ticks = [t[1],t[end]],
                y_ticks = [0,1],
                toggle_ticks_lab = [false,false],
               )
# Define vector of indices for plotting purposes
idx = collect(UnitRange(1::Int64, Nx))
randomised = collect(UnitRange(1::Int64, Nx))
StatsBase.direct_sample!(idx, randomised)
# Loop over the different ICs
for m in 1:Nx
        # Extract the solution
        x = X[:,m]
        # Plot the solution 
        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [(:teal,0.35), (:brown2,0.35)], colorrange = (1,Nx))
end
# Plot the equilibria
lines!(ax, [t[1], t[end]], [0.0, 0.0], linewidth = 6, color = :black)
lines!(ax, [t[1], t[end]], [1.0, 1.0], linewidth = 6, color = :black)
lines!(ax, [t[1], t[end]], [x3(A[1]) for T in [t[1], t[end]]], linewidth = 6, color = :black, linestyle = :dash)

#####################################
#   Dominant strategy (x1 stable)   #
#####################################

# Import the solution data at current timestep from csv 
solution = readin("../data/fig:autonomous/2.csv")
t = solution[:,1]
X = solution[:,2:end]

# Create and customise the figure 
fig, ax = mkfig(fig = fig,
                box_position = [1,1],
                limits = ((t[1], t[end]), (-0.05, 1.05)),
                title = L"\textbf{Dominant strategy (}\mathbf{x_1}\textbf{ stable)}",
                toggle_title = true,
                title_size = 30,
                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                toggle_lab = [false,true],
                lab_pad = [-60.0,-30.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [0,1],
                toggle_ticks_lab = [false,true],
                ticks_lab_xpos = [:right,:top],
                ticks_lab_trunc = [0,0],
               )
# Define vector of indices for plotting purposes
idx = collect(UnitRange(1::Int64, Nx))
randomised = collect(UnitRange(1::Int64, Nx))
StatsBase.direct_sample!(idx, randomised)
# Loop over the different ICs
for m in 1:Nx
        # Extract the solution
        x = X[:,m]
        # Plot the solution 
        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [(:teal,0.35), (:brown2,0.35)], colorrange = (1,Nx))
end
# Plot the equilibria
lines!(ax, [t[1], t[end]], [0.0, 0.0], linewidth = 6, color = :black)
lines!(ax, [t[1], t[end]], [1.0, 1.0], linewidth = 6, color = :black, linestyle = :dash)
lines!(ax, [t[1], t[end]], [x3(A[2]) for T in [t[1], t[end]]], linewidth = 6, color = :black, linestyle = :dash)

#####################################
#   Dominant strategy (x2 stable)   #
#####################################

# Import the solution data at current timestep from csv 
solution = readin("../data/fig:autonomous/3.csv")
t = solution[:,1]
X = solution[:,2:end]

# Create and customise the figure 
fig, ax = mkfig(fig = fig,
                box_position = [2,2],
                limits = ((t[1], t[end]), (-0.05, 1.05)),
                title = L"\textbf{Dominant strategy (}\mathbf{x_2}\textbf{ stable)}",
                toggle_title = true,
                title_size = 30,
                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                toggle_lab = [true,false],
                lab_pad = [-60.0,-30.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [0,1],
                toggle_ticks_lab = [true,false],
                ticks_lab_xpos = [:left,:top],
                ticks_lab_trunc = [0,0],
               )
# Define vector of indices for plotting purposes
idx = collect(UnitRange(1::Int64, Nx))
randomised = collect(UnitRange(1::Int64, Nx))
StatsBase.direct_sample!(idx, randomised)
# Loop over the different ICs
for m in 1:Nx
        # Extract the solution
        x = X[:,m]
        # Plot the solution 
        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [(:teal,0.35), (:brown2,0.35)], colorrange = (1,Nx))
end
# Plot the equilibria
lines!(ax, [t[1], t[end]], [0.0, 0.0], linewidth = 6, color = :black, linestyle = :dash)
lines!(ax, [t[1], t[end]], [1.0, 1.0], linewidth = 6, color = :black)
lines!(ax, [t[1], t[end]], [x3(A[3]) for T in [t[1], t[end]]], linewidth = 6, color = :black, linestyle = :dash)

#########################
#   Anti-coordination   #
#########################

# Import the solution data at current timestep from csv 
solution = readin("../data/fig:autonomous/4.csv")
t = solution[:,1]
X = solution[:,2:end]

# Create and customise the figure 
fig, ax = mkfig(fig = fig,
                box_position = [2,1],
                limits = ((t[1], t[end]), (-0.05, 1.05)),
                title = L"\textbf{Anti-coordination}",
                toggle_title = true,
                title_size = 30,
                lab = [L"\mathbf{t}",L"\mathbf{x}"],
                lab_pad = [-60.0,-30.0],
                x_ticks = [t[1],t[end]],
                y_ticks = [0,1],
                ticks_lab_xpos = [:right,:top],
                ticks_lab_trunc = [0,0],
               )
# Define vector of indices for plotting purposes
idx = collect(UnitRange(1::Int64, Nx))
randomised = collect(UnitRange(1::Int64, Nx))
StatsBase.direct_sample!(idx, randomised)
# Loop over the different ICs
for m in 1:Nx
        # Extract the solution
        x = X[:,m]
        # Plot the solution 
        lines!(ax, t, x, linewidth = 4, color = randomised[m], colormap = [(:teal,0.35), (:brown2,0.35)], colorrange = (1,Nx))
end
# Plot the equilibria
lines!(ax, [t[1], t[end]], [0.0, 0.0], linewidth = 6, color = :black, linestyle = :dash)
lines!(ax, [t[1], t[end]], [1.0, 1.0], linewidth = 6, color = :black, linestyle = :dash)
lines!(ax, [t[1], t[end]], [x3(A[4]) for T in [t[1], t[end]]], linewidth = 6, color = :black)

# Export the figure 
save("../fig/fig:autonomous_solutions.png", fig)
