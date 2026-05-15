include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# ==========================================
# | Discretisation settings
# ==========================================

Ω = [-1,1]    # Spatial domain 
BC = [1,2]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs 
Ωh = LinRange(Ω[1], Ω[2], Nh+1)   # Discretised domain
T = 1.0    # Time horizon 
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps

P1 = collect(LinRange(0.025, 0.1, 1000))   # Training set of problem 1
P2 = collect(LinRange(0.025, 0.1, 10))    # Training set of problem 2

# ==========================================
# | Problem 1
# ==========================================

# Import data from csv
X = readin("../data/heat_eq/problem_1/snapshots.csv")
S = readin("../data/heat_eq/problem_1/singular_values.csv")
V = readin("../data/heat_eq/problem_1/singular_vector.csv")
Vr = readin("../data/heat_eq/problem_1/pod_basis.csv")
coefficients = readin("../data/heat_eq/problem_1/coefficients.csv")
solutions = readin("../data/heat_eq/problem_1/solutions.csv")
u_fom = solutions[:,1]
u_rom = solutions[:,2]
r = size(Vr,2)

#################
#   Snapshots   #
#################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                box_position = [1:2,1],
                bg_out = :white,
                pad = (35,30,30,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (0.75,2.25)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-40.0,-20.0],
                lab_size = [30,30],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [1,2],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Loop over the columns of the snapshot matrix
for j in 1:size(X,2)
        # Reconstruct the FOM solution
        Uh = [BC[1]; X[:, j]; BC[2]]
        # Plot the snapshot
        lines!(ax, Ωh, Uh, linewidth = 5, color = j, colormap = (:thermal, 0.25), colorrange = (1,size(X,2)))
end
# Plot the colorbar
Colorbar(fig[1:2,2],
         size = 45, 
         spinewidth = 5.0,
         limits = (P1[1], P1[end]), 
         colormap = (:thermal, 0.5),
         ticks = [P1[1], P1[end]],
         ticksize = 22,
         tickwidth = 5.0,
         tickformat = "{:.$(2)f}",
         ticklabelsize = 30,
         label = L"\mathbf{\mu}",
         labelpadding = -40.0,
         labelsize = 30
        )

# Export the figure
save("../fig/heat_eq/fig:problem1_fom.png", fig)

#############################
#   Singular values decay   #
#############################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((0.75,5.25), (-100,1600)),
                lab = [L"\text{index}(j)", L"\mathbf{\sigma}_{j}"],
                lab_pad = [-60.0,-20.0],
                x_ticks = [1,5],
                y_ticks = [0,1550],
                ticks_lab_trunc = [0,0],
               )
# Plot the singular values decay
scatter!(ax, 1:10, S[1:10], markersize = 35, color = :red, strokewidth = 4, strokecolor = :black)
scatter!(ax, 1:r, S[1:r], markersize = 35, color = :green, strokewidth = 4, strokecolor = :black)

# Export the figure
save("../fig/heat_eq/fig:problem1_decay.png", fig)

##############################
#   Reduced order solution   #
##############################

# Create and customise the comparison figure 
fig, ax = mkfig(size = [2400,1600],
                box_position = [1:2,1:r],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (0.75,2.25)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [1,2],
                ticks_lab_trunc = [0,0],
               )
# Plot the FOM solution
lines!(ax, Ωh, u_fom, linewidth = 7, color = (:red,1))
# Plot the ROM solution
lines!(ax, Ωh, u_rom, linewidth = 10, linestyle = :dash, color = (:blue,1))

# Loop over the POD-basis vectors
for n in 1:r
        # Extract the coefficient of the ROM
        u_bar = trunc(coefficients[n], digits=2)
        # Create and customise the POD-basis figure 
        nullfig, ax = mkfig(fig=fig,
                            box_position = [3,n],
                            limits = ((Ω[1],Ω[end]), (-0.25,2.25)),
                            title = L"\bar{u}_{%$(n)} = %$(u_bar)",
                            toggle_title = true,
                            title_size = 50,
                            title_color = :black,
                            title_gap = 4.0,
                            lab = [L"\mathbf{x}", L"\mathbf{v}_{%$(n)}"],
                            lab_pad = [-60.0,-40.0],
                            x_ticks = [Ω[1],Ω[end]],
                            y_ticks = [0,2],
                            ticks_lab_trunc = [0,0],
                           )
        # Plot the basis function
        lines!(ax, Ωh[2:(end-1)], Vr[:,n], linewidth = 5, color = (:red,1))
        # Plot the weighted basis function
        lines!(ax, Ωh[2:(end-1)], u_bar.*Vr[:,n], linewidth = 7, linestyle = :dash, color = (:blue,1))
end

# Export the figure
save("../fig/heat_eq/fig:problem1_rom.png", fig)

# ==========================================
# | Problem 2
# ==========================================

# Import data from csv
X = readin("../data/heat_eq/problem_2/snapshots.csv")
S = readin("../data/heat_eq/problem_2/singular_values.csv")
V = readin("../data/heat_eq/problem_2/singular_vector.csv")
Vr = readin("../data/heat_eq/problem_2/pod_basis.csv")
coefficients = readin("../data/heat_eq/problem_2/coefficients.csv")
u_fom = readin("../data/heat_eq/problem_2/FOM_solutions.csv")
u_rom = readin("../data/heat_eq/problem_2/ROM_solutions.csv")
r = size(Vr,2)

#################
#   Snapshots   #
#################

# Create and customise the figure for μ1 
fig, ax = mkfig(size = [1200,800],
                box_position = [2,1],
                bg_out = :white,
                pad = (35,30,30,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (0.75,2.25)),
                title = L"\mathbf{\mu}_{1} = 0.1",
                toggle_title = true,
                title_size = 30,
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-40.0,-20.0],
                lab_size = [30,30],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [1,2],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Loop over the columns of the snapshot matrix
for j in 1:100
        # Reconstruct the FOM solution
        Uh = [BC[1]; X[:, j]; BC[2]]
        # Plot the snapshot
        lines!(ax, Ωh, Uh, linewidth = 5, color = j, colormap = (:thermal, 0.25), colorrange = (1,101))
end

# Create and customise the figure for μ10
fig, ax = mkfig(fig = fig,
                box_position = [2,2],
                limits = ((Ω[1],Ω[end]), (0.75,2.25)),
                title = L"\mathbf{\mu}_{10} = 0.25",
                toggle_title = true,
                title_size = 30,
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                toggle_lab = [true, false],
                lab_pad = [-40.0,-20.0],
                lab_size = [30,30],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [1,2],
                toggle_ticks_lab = [true, false],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,0],
               )
# Loop over the columns of the snapshot matrix
for j in 901:1001
        # Reconstruct the FOM solution
        Uh = [BC[1]; X[:, j]; BC[2]]
        # Plot the snapshot
        lines!(ax, Ωh, Uh, linewidth = 5, color = j, colormap = (:thermal, 0.25), colorrange = (901,1001))
end

# Plot the colorbar
Colorbar(fig[1,1:2],
         vertical = false,
         size = 45, 
         spinewidth = 5.0,
         limits = (P2[1], P2[end]), 
         colormap = (:thermal, 0.5),
         ticks = [P2[1], P2[end]],
         ticksize = 22,
         tickwidth = 5.0,
         tickformat = "{:.$(2)f}",
         ticklabelsize = 30,
         label = L"\mathbf{\mu}",
         labelpadding = -40.0,
         labelsize = 30
        )

# Export the figure
save("../fig/heat_eq/fig:problem2_fom.png", fig)

#############################
#   Singular values decay   #
#############################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((0.75,5.25), (-100,1600)),
                lab = [L"\text{index}(j)", L"\mathbf{\sigma}_{j}"],
                lab_pad = [-60.0,-20.0],
                x_ticks = [1,5],
                y_ticks = [0,1550],
                ticks_lab_trunc = [0,0],
               )
# Plot the singular values decay
scatter!(ax, 1:10, S[1:10], markersize = 35, color = :red, strokewidth = 4, strokecolor = :black)
scatter!(ax, 1:2, S[1:2], markersize = 35, color = :green, strokewidth = 4, strokecolor = :black)

# Export the figure
save("../fig/heat_eq/fig:problem2_decay.png", fig)

##############################
#   Reduced order solution   #
##############################

# Define timesteps for plotting
T_plt = [1,10,100,1000]

# Create and customise the comparison figure 
fig, ax = mkfig(size = [2400,1600],
                box_position = [1:2,1:2],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((Ω[1],Ω[end]), (0.75,2.25)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [1,2],
                ticks_lab_trunc = [0,0],
               )
# Loop over the timesteps
for t in T_plt
        # Plot the FOM solution
        lines!(ax, Ωh, u_fom[:,t], linewidth = 7, color = (:red,1))
        # Plot the ROM solution
        lines!(ax, Ωh, u_rom[:,t], linewidth = 10, linestyle = :dash, color = (:blue,1))
end

# Extract the reduce system's  solution
u_bar_1 = coefficients[1,:]
u_bar_2 = coefficients[2,:]

# Create and customise the 1st POD-basis figure 
fig, ax = mkfig(fig=fig,
                box_position = [3,1],
                limits = ((Ω[1],Ω[end]), (0.02,0.042)),
                lab = [L"\mathbf{x}", L"\mathbf{v}_{1}"],
                toggle_lab = [false,true],
                toggle_ticks_lab = [false,true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [0.02,0.042],
                ticks_lab_trunc = [0,2],
               )
# Plot the basis function
lines!(ax, Ωh[2:(end-1)], Vr[:,1], linewidth = 5, color = (:red,1))

# Create and customise the 1st reduced coefficients solution figure 
fig, ax = mkfig(fig=fig,
                box_position = [3,2],
                limits = ((0,T), (48,50)),
                lab = [L"\mathbf{t}", L"\mathbf{\bar{u}}_1"],
                toggle_lab = [false,true],
                toggle_ticks_lab = [false,true],
                lab_pad = [-60.0,-40.0],
                x_ticks = [0,T],
                y_ticks = [48,50],
                ticks_lab_trunc = [0,0],
               )
# Plot the reduced coefficients time evolution 
lines!(ax, LinRange(0, T, Nt+1), u_bar_1, linewidth = 5, linestyle = :dash, color = (:blue,1))

# Create and customise the 2nd POD-basis figure 
fig, ax = mkfig(fig=fig,
                box_position = [4,1],
                limits = ((Ω[1],Ω[end]), (-0.045,0.05)),
                lab = [L"\mathbf{x}", L"\mathbf{v}_{2}"],
                lab_pad = [-60.0,-75.0],
                x_ticks = [Ω[1],Ω[end]],
                y_ticks = [-0.04,0.05],
                ticks_lab_trunc = [0,2],
               )
# Plot the basis function
lines!(ax, Ωh[2:(end-1)], Vr[:,2], linewidth = 5, color = (:red,1))

# Create and customise the 1st reduced coefficients solution figure 
fig, ax = mkfig(fig=fig,
                box_position = [4,2],
                limits = ((0,T), (-20,5)),
                lab = [L"\mathbf{t}", L"\mathbf{\bar{u}}_2"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [0,T],
                y_ticks = [-20,5],
                ticks_lab_trunc = [0,0],
               )
# Plot the reduced coefficients time evolution 
lines!(ax, LinRange(0, T, Nt+1), u_bar_2, linewidth = 5, linestyle = :dash, color = (:blue,1))



# Export the figure
save("../fig/heat_eq/fig:problem2_rom.png", fig)
