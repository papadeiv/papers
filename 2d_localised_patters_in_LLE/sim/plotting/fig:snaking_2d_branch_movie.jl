include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

using MAT

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

# Extract the indices of the distinct branches 
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

# Define current branch
lower_idx = (branch_idx_plt[2]÷10)*10
upper_idx = (branch_idx_plt[3]÷10)*10
current_branch = collect(lower_idx:10:upper_idx)
N_steps = length(current_branch) 

# Define a vector of transparency values for the branches
branch_alpha = [0.25,1.0,0.25]

# Define variables and parameters
L = 80.0::Float64
x = collect(LinRange(-L,L,1000)) 
Nh = 500::Int64

# Create box 1 in the main bifurcation diagram
ρ_1_min = 1.1175 
ρ_1_max = 1.1180 
u_1_min = 0.0250
u_1_max = 0.0650
box_1 = [Point(ρ_1_min,u_1_min),
         Point(ρ_1_max,u_1_min),
         Point(ρ_1_max,u_1_max),
         Point(ρ_1_min,u_1_max),
         Point(ρ_1_min,u_1_min),
        ]

# Create box 2 in the main bifurcation diagram
ρ_2_min = 1.1150
ρ_2_max = 1.1170
u_2_min = 0.0500
u_2_max = 0.0750
box_2 = [Point(ρ_2_min,u_2_min),
         Point(ρ_2_max,u_2_min),
         Point(ρ_2_max,u_2_max),
         Point(ρ_2_min,u_2_max),
         Point(ρ_2_min,u_2_min),
        ]

# Loop over continuation steps of the branch
#=@showprogress=#for n in 1:N_steps
        # Create and customise the main bifurcation diagram subplot 
        fig, bd_ax = mkfig(size = [1200,1200],
                           bg_out = :white,
                           box_position = [1:2,1:2],
                           pad = (5,50,25,25), # Order is: left, right, bottom, top 
                           limits = ((1.111, 1.119), (0, 0.17)),
                           lab = [L"\mathbf{\rho}", L"\mathbf{||u - u_0||}_{2}"],
                           lab_pad = [-30.0,-30.0],
                           lab_size = [30,30],
                           x_ticks = [1.111,1.119],
                           y_ticks = [minimum(u_2),maximum(u_2)],
                           ticks_lab_size = [30,30],
                           ticks_lab_trunc = [3,2],
                          )
        # Loop over the branches
        for n in 1:N_branches
                # Plot the branch
                lines!(bd_ax, branches[n][:,1], branches[n][:,2], linewidth = 3.0, color = (:magenta4, branch_alpha[n]))
        end
        # Plot the branch of the second continuation
        lines!(bd_ax, ρ_2, u_2_2, linewidth = 3.0, color = (:magenta4, 0.25))
        # Plot box 1 for the 1st zoomed-in subplot and its text label
        poly!(bd_ax, box_1, color=:transparent, strokecolor=:black, strokewidth=2)
        text!(ρ_1_max, u_1_max, text = L"\mathbf{A}", align = (:left, :bottom), fontsize = 30)
        # Plot box 2 for the 2nd zoomed-in subplot and its text label
        poly!(bd_ax, box_2, color=:transparent, strokecolor=:black, strokewidth=2)
        text!(ρ_2_max, u_2_max, text = L"\mathbf{B}", align = (:left, :bottom), fontsize = 30, color = :black)

        # Create and customise the 1st zoomed-in subplot (box 1) 
        nullfig, bd_ax_1 = mkfig(fig = fig,
                                 box_position = [3:4,1:2],
                                 limits = ((ρ_1_min,ρ_1_max), (u_1_min,u_1_max)),
                                 toggle_lab = [false,false],
                                 lab_pad = [0.0,0.0],
                                 lab_size = [30,30],
                                 x_ticks = [ρ_1_min,ρ_1_max],
                                 y_ticks = [u_1_min,u_1_max],
                                 ticks_lab_size = [30,30],
                                 ticks_lab_trunc = [4,3],
                                )
        # Loop over the branches
        for n in 1:N_branches
                # Plot the branch
                lines!(bd_ax_1, branches[n][:,1], branches[n][:,2], linewidth = 3.0, color = (:magenta4, branch_alpha[n]))
        end
        # Plot the text label
        text!(1.11795, 0.0625, text = L"\mathbf{A}", align = (:right, :top), fontsize = 30, color = :black)

        # Create and customise the 2nd zoomed-in subplot 
        nullfig, bd_ax_2 = mkfig(fig = fig,
                                 box_position = [1:2,3:4],
                                 limits = ((ρ_2_min,ρ_2_max), (u_2_min,u_2_max)),
                                 toggle_lab = [false,false],
                                 lab_pad = [0.0,0.0],
                                 lab_size = [30,30],
                                 x_ticks = [ρ_2_min,ρ_2_max],
                                 y_ticks = [u_2_min,u_2_max],
                                 ticks_lab_size = [30,30],
                                 ticks_lab_trunc = [4,3],
                                )
        # Loop over the branches
        for n in 1:N_branches
                # Plot the branch
                lines!(bd_ax_2, branches[n][:,1], branches[n][:,2], linewidth = 3.0, color = (:magenta4, branch_alpha[n]))
        end
        # Plot the text label
        text!(1.11685, 0.0735, text = L"\mathbf{B}", align = (:right, :top), fontsize = 30, color = :black)

        display(n)
        display(current_branch[n])

        # Extract the continuation parameter value and the L2-norm
        ρ_current = (branches[2])[current_branch[n],1]
        u_2_current = (branches[2])[current_branch[n],2]

        # Plot the solution's norm in the various bifurcation diagrams
        scatter!(bd_ax, ρ_current, u_2_current, markersize = 30, strokewidth = 5, strokecolor = :black, color = :lightslateblue)
        scatter!(bd_ax_1, ρ_current, u_2_current, markersize = 30, strokewidth = 5, strokecolor = :black, color = :lightslateblue)
        scatter!(bd_ax_2, ρ_current, u_2_current, markersize = 30, strokewidth = 5, strokecolor = :black, color = :lightslateblue)

        # Read in solution from formatted index
        filename = @sprintf("../data/cont_rho_2d/solution_%07d.mat", current_branch[n])
        file = matread(filename)
        
        # Extract the real solution field
        u = file["u"]
        u_re = u[1:Nh]

        # Define a discretised polar grid
        r = range(0, L, length=Nh) 
        θ = range(0, 2*pi, length=Nh)

        # Reconstruct the 2-dimensional solution in polar coordinates
        u_re_2d = [u_re[k] for j in θ, k in 1:length(r)]

        # Create the polar grid subplot
        ax = PolarAxis(fig[3:4,3:4],
                       rlimits = (:origin, L),
                       rgridvisible = false,
                       rticklabelsvisible = false,
                       thetagridvisible = false,
                       thetaticklabelsvisible = false,
                       spinewidth = 5.0,
                      )
        # Plot the contour of the 2-dimensional field
        surface!(ax, 0..2pi, 0..L, u_re_2d, color = u_re_2d, shading = NoShading, colormap = :viridis)
        tightlimits!(ax)

        # Adjust whitespace between rows and columns
        rowgap!(fig.layout, 50)

        # Export the figure 
        save("../fig/fig:snaking_2d/$n.png", fig)
end
