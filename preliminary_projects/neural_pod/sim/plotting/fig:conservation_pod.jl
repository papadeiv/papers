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
X = readin("../data/conservation_law/snapshots.csv")
S = readin("../data/conservation_law/singular_values.csv")
V = readin("../data/conservation_law/singular_vector.csv")
Vr = readin("../data/conservation_law/pod_basis.csv")
coefficients = readin("../data/conservation_law/coefficients.csv")
u_fom = readin("../data/conservation_law/FOM_solutions.csv")
u_rom = readin("../data/conservation_law/ROM_solutions.csv")
r = size(Vr,2)

#################
#   Snapshots   #
#################

# Define timesteps for plotting
T_plt = [1,50,100,150,200,400,600]

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
save("../fig/fig:conservation_fom.png", fig)

#############################
#   Singular values decay   #
#############################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((0.75,10.25), (-10,180)),
                lab = [L"\text{index}(j)", L"\mathbf{\sigma}_{j}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [1,10],
                y_ticks = [0,160],
                ticks_lab_trunc = [0,0],
               )
# Plot the singular values decay
scatter!(ax, 1:10, S[1:10], markersize = 35, color = :red, strokewidth = 4, strokecolor = :black)
scatter!(ax, 1:r, S[1:r], markersize = 35, color = :green, strokewidth = 4, strokecolor = :black)

# Export the figure
save("../fig/fig:conservation_decay.png", fig)

##############################
#   Reduced order solution   #
##############################

# Define timesteps for plotting
T_plt = [1,200,400,600]

# Create and customise the comparison figure 
fig, ax = mkfig(size = [2400,1600],
                box_position = [1:2,1:r],
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
        lines!(ax, Ωh, [BC[1]; u_fom[:,t]; BC[2]], linewidth = 7, color = (:red,1))
        # Plot the ROM solution
        lines!(ax, Ωh, u_rom[:,t], linewidth = 10, linestyle = :dash, color = (:blue,1))
end

# Loop over the POD-basis vectors
for n in 1:r
        # Create and customise the POD-basis figure 
        nullfig, ax = mkfig(fig=fig,
                            box_position = [3,n],
                            limits = ((Ω[1],Ω[end]), (-0.125,0.125)),
                            lab = [L"\mathbf{x}", L"\mathbf{v}_{%$(n)}"],
                            lab_pad = [-60.0,-60.0],
                            x_ticks = [Ω[1],Ω[end]],
                            y_ticks = [-0.1,0.1],
                            ticks_lab_trunc = [0,1],
                           )
        # Plot the basis function
        lines!(ax, Ωh[2:(end-1)], Vr[:,n], linewidth = 7, color = n, colormap = :darktest, colorrange = (1,r))
end

# Create and customise the reduced coefficients figure 
fig, ax = mkfig(fig=fig,
                box_position = [4:5,1:r],
                limits = ((0,T), (-6,7)),
                lab = [L"\mathbf{t}", L"\bar{\mathbf{u}}_{j}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [0,T],
                y_ticks = [-6,7],
                ticks_lab_trunc = [0,0],
               )
# Loop over the reduced coefficients 
for n in 1:r
        # Plot the reduced coefficients 
        lines!(ax, LinRange(0, T, Nt+1), coefficients[n,:], linewidth = 7, color = n , colormap = :darktest, colorrange = (1,r))
end

# Export the figure
save("../fig/fig:conservation_rom.png", fig)
