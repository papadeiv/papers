include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

using MAT

################################
#     Bifurcation diagram      #
################################

# Import data from mat files
file = matread("../data/cont_rho_2d/branch.mat")
key = file["branch"]
λ = key[:,2]
ρ = key[:,3]
u_2 = key[:,5]

file = matread("../data/cont_rho_2d_2/branch.mat")
key = file["branch"]
λ_2 = key[:,2]
ρ_2 = key[:,3]
u_2_2 = key[:,5]

# Define plotting indices
contiguous_zeros = findall(x -> isapprox(x, 0.0, atol=1e-3), u_2)
branch_idx_plt = contiguous_zeros[[true; diff(contiguous_zeros) .!= 1]]
branch_idx_plt = [1; branch_idx_plt ; length(ρ)]
N_branches = length(branch_idx_plt)-1

# Extract the individual branches
branches = Vector{Matrix{Float64}}(undef, N_branches)
for n in 1:N_branches
        # Get the number of steps in the current branch
        N_steps = (branch_idx_plt[n+1] - branch_idx_plt[n]) + 1
        # Define a placehold matrix
        A = Matrix{Float64}(undef, N_steps, 2)
        # Extract the parameter values
        A[:,1] = ρ[branch_idx_plt[n]:branch_idx_plt[n+1]]
        # Extract the norm values
        A[:,2] = u_2[branch_idx_plt[n]:branch_idx_plt[n+1]]
        # Fill-in the branch
        branches[n] = A 
end

# Define a vector of transparency values for the branches
branch_alpha = [0.25,1.0,0.25,0.25]

# Define number of columns for the functions subplots
N_plt = 3::Int64

#=
# Extract the indices of the saddle-node bifurcations 
idx_bif = findall(diff(λ) .!=0)

# Check whether there have been bifurcations
if length(idx_bif) > 0
        # Separate the stable solutions from the unstable ones
        global idx_stb = Vector{Float64}(undef, length(idx_bif))

        # Update the stability array
        for n in 1:length(idx_bif)
                if λ[idx_bif[n]-1] > 0
                        idx_stb[n] = +1
                else
                        idx_stb[n] = -1
                end
        end

        # Update the bifurcation and stability indices to include the endpoints
        idx_bif = [1; idx_bif; length(λ)]
        idx_stb = [idx_stb[2]; idx_stb; idx_stb[end]]

        # Define a global variable for the number of bifurcations
        global N_bif = length(idx_bif)-1
else
        # Create a bifurcation index vector that spans the whole branch  
        global idx_bif = [1, length(λ)]
        
        # Create a single-entry stability vector
        if λ[1] > 0
                global idx_stb = [+1]
        else
                global idx_stb = [-1]
        end

        # Define a global variable for the number of bifurcations
        global N_bif = 1
end
=#

# Create and customise the figure 
fig, ax = mkfig(size = [1200,1200],
                bg_out = :white,
                box_position = [1:3,1:N_plt],
                pad = (5,30,25,25), # Order is: left, right, bottom, top 
                limits = ((1.110, 1.119), (-0.005, 0.275)),
                lab = [L"\mathbf{\rho}", L"\mathbf{||u - u_0||}_{2}"],
                lab_pad = [0.0,0.0],
                lab_size = [30,30],
                x_ticks = [1.110,1.11217,1.11452,1.119],
                y_ticks = [0,0.27],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [3,2],
               )
# Plot vertical lines for the folding points in 1 dimension
lines!(ax, [1.11208,1.11208], [-0.005,0.275], linewidth = 1.0, color = :black, linestyle = :dash)
lines!(ax, [1.11467,1.11467], [-0.005,0.275], linewidth = 1.0, color = :black, linestyle = :dash)
# Loop over the branches
for n in 1:N_branches
        # Plot the branch
        lines!(ax, branches[n][:,1], branches[n][:,2], linewidth = 3.0, color = (:magenta4, branch_alpha[n]))
end
# Plot the branch of the second continuation run
lines!(ax, ρ_2, u_2_2, linewidth = 3.0, color = (:magenta4, 0.25))
# Plot vertical lines for the folding points in 2 dimensions
lines!(ax, [1.11217,1.11217], [-0.005,0.275], linewidth = 1.0, color = :black)
lines!(ax, [1.11452,1.11452], [-0.005,0.275], linewidth = 1.0, color = :black)

# Loop over the saddle-node bifurcation indices
#=
for n in 1:(length(idx_bif)-1)
        # Check the stability of the subbranch
        if idx_stb[n] < 0
                # Plot the branch of stable solutions
                lines!(ax, ρ[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 3.0, color = :magenta4)
        else
                # Plot the branch of unstable solutions
                lines!(ax, ρ[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 3.0, color = :magenta4, linestyle = :dash)
        end
end
=#

#####################
#     Solutions     #
#####################

# Define variables and parameters
L = 80.0::Float64
x = collect(LinRange(-L,L,1000)) 
Nh = 500::Int64
U = Matrix{Float64}(undef, Nh, 2*N_plt)
ρ_plt = Vector{Float64}(undef, 2*N_plt)
idx_plt = [1,61,301,1001,2101,2901] 

# Loop over the plotting indices
for n in 1:(2*N_plt)
        # Read in solution from formatted index
        filename = @sprintf("../data/cont_rho_2d/solution_%07d.mat", idx_plt[n]-1)
        file = matread(filename)
        
        # Extract the real solution field
        u = file["u"]
        U[:,n] = u[1:Nh]

        # Extract the continuation parameter value
        ρ = file["p"]
        ρ_plt[n] = ρ[2]
end

# Locate the solutions you're about to plot in the bifurcation diagram
scatter!(ax, ρ_plt, u_2[idx_plt], markersize = 30, strokewidth = 5, strokecolor = :black, color = :lightslateblue)
# Locate the IC of the continuation problem 
scatter!(ax, ρ_plt[1], u_2[idx_plt[1]], marker = :star5, markersize = 40, strokewidth = 3, strokecolor = :black, color = :yellow)

# Define the size of the polar subplots
plot_size = 200

# First row of solutions
for n in 1:N_plt
        # Extract the solution field
        u_re = U[:,n]

        # Define a discretised polar grid
        r = range(0, L, length=Nh) 
        θ = range(0, 2*pi, length=Nh)

        # Reconstruct the 2-dimensional solution in polar coordinates
        u_re_2d = [u_re[m] for n in θ, m in 1:length(r)]

        # Create the polar grid subplot
        ax = PolarAxis(fig[4,n], 
                       rlimits = (:origin, L),
                       rgridvisible = false,
                       rticklabelsvisible = false,
                       thetagridvisible = false,
                       thetaticklabelsvisible = false,
                       spinewidth = 5.0,
                       width = plot_size,
                       height = plot_size,
                      )
        # Plot the contour of the 2-dimensional field
        global s = surface!(ax, 0..2pi, 0..L, u_re_2d, color = u_re_2d, shading = NoShading, colormap = :viridis)
        colsize!(fig.layout, n, Auto(1)) 
        tightlimits!(ax)
end

# Second row of solutions
for n in 1:N_plt
        # Extract the solution field
        u_re = U[:,n+N_plt]

        # Define a discretised polar grid
        r = range(0, L, length=Nh) 
        θ = range(0, 2*pi, length=Nh)

        # Reconstruct the 2-dimensional solution in polar coordinates
        u_re_2d = [u_re[m] for n in θ, m in 1:length(r)]

        # Create the polar grid subplot
        ax = PolarAxis(fig[5,n], 
                       rlimits = (:origin, L),
                       rgridvisible = false,
                       rticklabelsvisible = false,
                       thetagridvisible = false,
                       thetaticklabelsvisible = false,
                       spinewidth = 5.0,
                       width = plot_size,
                       height = plot_size,
                      )
        # Plot the contour of the 2-dimensional field
        global s = surface!(ax, 0..2pi, 0..L, u_re_2d, color = u_re_2d, shading = NoShading, colormap = :viridis)
        colsize!(fig.layout, n, Auto(1)) 
        tightlimits!(ax)
end

# Plot the colorbar
Colorbar(fig[6,1:N_plt],
         size = 45,
         spinewidth = 5.0,
         vertical = false,
         ticks = [minimum(U),maximum(U)],
         limits = [minimum(U),maximum(U)],
         label = L"\textbf{Re}\mathbf{(u)}",
         labelsize = 30,
         labelpadding = -40.0,
         ticklabelsize = 30,
         ticksize = 22, 
         tickwidth = 5.0,
         tickformat = "{:.1f}" 
        )

# Adjust whitespace between rows and columns
colgap!(fig.layout, 200)
rowgap!(fig.layout, 30)

# Export the figure 
save("../fig/fig:snaking_2d.png", fig)
