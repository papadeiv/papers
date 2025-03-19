using MAT

# Creat the figure for the BD
using CairoMakie, Makie.Colors
CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
ax = Axis(fig[1, 1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-1, 1.15), (-1, 20)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 25,
    xlabelcolor = :black,
    xlabelpadding = 0,
    xticks = [-1,-0.9,0,1,1.1],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 25,
    xtickformat = "{:.1f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"||u||_{2}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = -30.0,
    yticks = [0,20],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = true,
    yticklabelsize = 25,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)

###############################################
#                GROUND STATE                 #
###############################################

# Print the variables' names
file = matopen("pt625.mat")
names = keys(file)
close(file)

# Extract the relevant keys
variables = matread("pt625.mat")
key = variables["p"]

# Extract the branch information
branch = key["branch"]
λ = branch[3,:]
μ = branch[4,:]
u = branch[5,:]

# Divide the branches by stability
μ_stable = Vector{Float64}()
u_stable = Vector{Float64}()
μ_unstable = Vector{Float64}()
u_unstable = Vector{Float64}()

for n in 1:length(λ)
        if λ[n] > 0.0
                push!(μ_unstable, μ[n])
                push!(u_unstable, u[n])
        else
                push!(μ_stable, μ[n])
                push!(u_stable, u[n])
        end
end

# Plot the branches
lines!(ax, μ_stable, u_stable, color = :black, linewidth = 2)
lines!(ax, μ_unstable, u_unstable, color = :black, linewidth = 2, linestyle = :dash)
# Store the data for the animation later
μ1 = μ_stable
u1 = u_stable
μ2 = μ_unstable
u2 = u_unstable

###############################################
#           HOMOGENEOUS SOLUTIONS             #
###############################################

# Extract the relevant keys
variables = matread("pt270.mat")
key = variables["p"]

# Extract the branch information of the first solution
branch = key["branch"]
λ = branch[3,:]
μ = branch[4,:]
u = branch[6,:]

# Divide the branches by stability
μ_stable = Vector{Float64}()
u_stable = Vector{Float64}()
μ_unstable = Vector{Float64}()
u_unstable = Vector{Float64}()

for n in 1:length(λ)
        if λ[n] > 0.0
                push!(μ_unstable, μ[n])
                push!(u_unstable, u[n])
        else
                push!(μ_stable, μ[n])
                push!(u_stable, u[n])
        end
end

# Plot the branches
lines!(ax, μ_stable, u_stable, color = :black, linewidth = 2)
lines!(ax, μ_unstable, u_unstable, color = :black, linewidth = 2, linestyle = :dash)
# Store the data for the animation later
μ3 = μ_stable
u3 = u_stable
μ4 = μ_unstable
u4 = u_unstable

###############################################
#         SPATIALLY-PERIODIC SOLUTION         #
###############################################

# Extract the relevant keys
variables = matread("pt393.mat")
key = variables["p"]

# Extract the branch information of the first solution
branch = key["branch"]
λ = branch[3,2:end]
μ = branch[4,2:end]
u = branch[6,2:end]

# Divide the branches by stability
μ_stable = Vector{Float64}()
u_stable = Vector{Float64}()
μ_unstable = Vector{Float64}()
u_unstable = Vector{Float64}()

for n in 1:length(λ)
        if λ[n] > 0.0
                push!(μ_unstable, μ[n])
                push!(u_unstable, u[n])
        else
                push!(μ_stable, μ[n])
                push!(u_stable, u[n])
        end
end

# Plot the branches
lines!(ax, μ_stable, u_stable, color = :blue, linewidth = 2)
lines!(ax, μ_unstable, u_unstable, color = :blue, linewidth = 2, linestyle = :dash)
# Store the data for the animation later
μ5 = μ_stable
u5 = u_stable
μ6 = μ_unstable
u6 = u_unstable

###############################################
#             TRANSIENT SOLUTION              #
###############################################

# Read the solution 
using CSV, DataFrames
df = DataFrame(CSV.File("./solution.csv"; delim=','))
μ_h = df[!,1]
u_2 = df[!,2]

# Plot the transient solution's norm 
#L = 64.0*π
L = 2.235
lines!(ax, μ_h, (1/L).*u_2, color = (:red, 0.25), linewidth = 3)

# Plot points on the transient
scatter!(ax, μ_h[1], (1/L)*u_2[1], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 17)
#scatter!(ax, μ_h[2075], u_2[2075], color = :yellow, strokecolor = :black, strokewidth = 1.5, markersize = 17)
scatter!(ax, μ_h[end], (1/L)*u_2[end], color = :green, strokecolor = :black, strokewidth = 1.5, markersize = 17)

# Export the BD
save("../../../results/pdes/SH35FS/fig3.3.png", fig)

# Create a flipbook animation
for n in 1:length(μ_h)
        local fig = Figure(; size = (1200, 800), backgroundcolor = "#fdfff2ff")
        local ax = Axis(fig[1, 1],
                # Background
                backgroundcolor = "#fdfff2ff",
                xgridvisible = false,
                ygridvisible = false,
                limits = ((-1, 1.15), (-1, 20)),
                # Title
                title = "Title",
                titlevisible = false,
                titlesize = 25,
                titlealign = :center,
                titlegap = -38.0,
                # x-axis
                xlabel = L"\mu",
                xlabelvisible = true,
                xlabelsize = 25,
                xlabelcolor = :black,
                xlabelpadding = 0,
                xticks = [-1,-0.9,0,1,1.1],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 25,
                xtickformat = "{:.1f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"||u||_{2}",
                ylabelvisible = true,
                ylabelsize = 30,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                yticks = [0,20],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 25,
                ytickformat = "{:.0f}",
                yscale = identity,
                yaxisposition = :left,
        )
        lines!(ax, μ1, u1, color = :black, linewidth = 2)
        lines!(ax, μ2, u2, color = :black, linewidth = 2, linestyle = :dash)
        lines!(ax, μ3, u3, color = :black, linewidth = 2)
        lines!(ax, μ4, u4, color = :black, linewidth = 2, linestyle = :dash)
        lines!(ax, μ5, u5, color = :blue, linewidth = 2)
        lines!(ax, μ6, u6, color = :blue, linewidth = 2, linestyle = :dash)
        lines!(ax, μ_h[1:n], (1/L).*u_2[1:n], color = (:red, 0.25), linewidth = 3)
        scatter!(ax, μ_h[n], (1/L)*u_2[n], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 17)
        save("../../../results/pdes/SH35FS/flipbook/bd/$n.png", fig)
end
