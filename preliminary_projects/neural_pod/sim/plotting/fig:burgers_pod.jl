include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# ==========================================
# | Discretisation settings
# ==========================================

Ω = [-1,1]    # Spatial domain 
BC = [0,0]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs 
Ωh = LinRange(Ω[1], Ω[2], Nh+1)   # Discretised domain
T = 1.0    # Time horizon 
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps

# Import data from csv
X = readin("../data/burgers_eq/snapshots.csv")
S = readin("../data/burgers_eq/singular_values.csv")
V = readin("../data/burgers_eq/singular_vector.csv")
Vr = readin("../data/burgers_eq/pod_basis.csv")
coefficients = readin("../data/burgers_eq/coefficients.csv")
u_fom = readin("../data/burgers_eq/FOM_solutions.csv")
u_rom = readin("../data/burgers_eq/ROM_solutions.csv")
r = size(Vr,2)

println("--- Generating the figures ---")

#################
#   Snapshots   #
#################

# Define timesteps for plotting
T_plt = [1,50,100,150,200,300,400,500,600,700,800,900,1000]

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (35,30,30,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (-0.15,1.15)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-40.0,-40.0],
                lab_size = [30,30],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [0,1],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Loop over the plotting indices 
for j in 1:length(T_plt)
        # Reconstruct the FOM solution
        Uh = [BC[1]; X[:,T_plt[j]]; BC[2]]
        # Compute the alpha value (transparency)
        alpha = convert(Float64, j/length(T_plt))
        # Plot the snapshot
        lines!(ax, Ωh, Uh, linewidth = 5, color = (:brown2, alpha))
end

# Export the figure
save("../fig/fig:burgers_fom.png", fig)

#############################
#   Singular values decay   #
#############################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((0,31), (-10,180)),
                lab = [L"\text{index}(j)", L"\mathbf{\sigma}_{j}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [1,30],
                y_ticks = [0,160],
                ticks_lab_trunc = [0,0],
               )
# Plot the singular values decay
scatter!(ax, 1:30, S[1:30], markersize = 35, color = :red, strokewidth = 4, strokecolor = :black)
scatter!(ax, 1:r, S[1:r], markersize = 35, color = :green, strokewidth = 4, strokecolor = :black)

# Export the figure
save("../fig/fig:burgers_decay.png", fig)

##############################
#   Reduced order solution   #
##############################

# Define timesteps for plotting
T_plt = [1,500,1000]

# Create and customise the comparison figure 
fig, ax = mkfig(size = [2400,1600],
                box_position = [1:2,1:5],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (-0.15,1.15)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [0,1],
                ticks_lab_trunc = [0,0],
               )
# Loop over the timesteps
for t in T_plt
        # Plot the FOM solution
        lines!(ax, Ωh, [BC[1]; u_fom[:,t]; BC[2]], linewidth = 10, color = (:red,1))
        # Plot the ROM solution
        lines!(ax, Ωh, u_rom[:,t], linewidth = 6, color = (:blue,1))
end

# Loop over the POD-basis vectors (first row)
for n in 1:5
        # Create and customise the POD-basis figure 
        local nullfig, ax = mkfig(fig=fig,
                            box_position = [3,n],
                            limits = ((Ω[1],Ω[end]), (-0.125,0.125)),
                            lab = [L"\mathbf{x}", L"\mathbf{v}_{%$(n)}"],
                            lab_pad = [-60.0,-60.0],
                            x_ticks = [Ω[1],Ω[end]],
                            y_ticks = [-0.1,0.1],
                            ticks_lab_trunc = [0,1],
                           )
        # Plot the basis function
        lines!(ax, Ωh[2:(end-1)], Vr[:,n], linewidth = 7, color = n, colormap = :darktest, colorrange = (1,10))
end

# Loop over the POD-basis vectors (second row)
for m in 1:5
        # Define shifted index
        n = m + 6
        # Create and customise the POD-basis figure 
        local nullfig, ax = mkfig(fig=fig,
                            box_position = [4,m],
                            limits = ((Ω[1],Ω[end]), (-0.125,0.125)),
                            lab = [L"\mathbf{x}", L"\mathbf{v}_{%$(n)}"],
                            lab_pad = [-60.0,-60.0],
                            x_ticks = [Ω[1],Ω[end]],
                            y_ticks = [-0.1,0.1],
                            ticks_lab_trunc = [0,1],
                           )
        # Plot the basis function
        lines!(ax, Ωh[2:(end-1)], Vr[:,n], linewidth = 7, color = n, colormap = :darktest, colorrange = (1,10))
end

# Export the figure
save("../fig/fig:burgers_rom.png", fig)
