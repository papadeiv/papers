include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

#############
# B-tipping #
#############

# Import the data from csv
u_B = readin("../data/slide_04/B-tipping.csv")
t = u_B[:,1] 
μ = u_B[:,2]
u = u_B[:,3]
stable_1 = readin("../data/slide_04/stable_eq_1.csv")
μ1 = stable_1[:,1]
x1 = stable_1[:,2]
stable_2 = readin("../data/slide_04/stable_eq_2.csv")
μ2 = stable_2[:,1]
x2 = stable_2[:,2]
unstable = readin("../data/slide_04/unstable_eq.csv")
μ3 = unstable[:,1]
x3 = unstable[:,2]

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

# Define the limits for the scalar potential plot
x_inf = -(1.0::Float64)
x_sup = 3.8::Float64
y_inf = -(2.0::Float64)
y_sup = 3.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Create and customise the explanatory potential function on the left
fig, ax = mkfig(size = [2000,1000],
                box_position = [1,1],
                border_color = RGBAf(0,0.236,0.236),
                pad = 30, # Equal amount of "whitespace" in all four directions 
                limits = ((x_inf, x_sup), (y_inf, y_sup)),
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ[302]) for x in domain], color = :black, linewidth = 4.5)
# Plot the state (ball in the landscape)
scatter!(ax, u[302], V(u[302],μ[302]), color = (RGBAf(0,0.236,0.236), 1.0), markersize = 45, strokewidth = 3)
 
# Create and customise the explanatory potential function on the right 
fig, ax = mkfig(fig = fig,
                box_position = [1,2],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((x_inf, x_sup), (y_inf, y_sup)),
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ[463]) for x in domain] .- 1.5, color = :black, linewidth = 4.5)
# Plot the state (ball in the landscape)
scatter!(ax, u[463], V(u[463],μ[463]) - 1.5, color = (RGBAf(0,0.236,0.236), 1.0), markersize = 45, strokewidth = 3)
 
# Export the phase plane figure 
save("../fig/slide_04_1.png", fig)

# Define the limits for the scalar potential plot
x_inf = -(3.0::Float64)
x_sup = 5.0::Float64
y_inf = -(12.5::Float64)
y_sup = 25.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Loop over the parameter values
printstyled("Generating figures for B-tipping using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:length(μ)
        # Create and customise the potential figure
        fig, ax = mkfig(size = [800,1200],
                        bg_out = colorant"rgb(147,237,213)",
                        border_color = RGBAf(0,0.236,0.236),
                        box_position = [2:3,1],
                        pad = 30, # Equal amount of "whitespace" in all four directions 
                        limits = ((x_inf, x_sup), (y_inf, y_sup)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[n]) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ[n]), color = (:teal, 1.0), markersize = 45, strokewidth = 3)
 
        # Create and customise the bifurcation diagram figure
        fig, ax = mkfig(fig = fig,
                        bg_out = colorant"rgb(147,237,213)",
                        border_color = RGBAf(0,0.236,0.236),
                        box_position = [1,1],
                        limits = ((μ[1], μ[end]), (-1, 3.5)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the bifurcation diagram 
        lines!(ax, μ1, x1, color = :black, linewidth = 3)
        lines!(ax, μ2, x2, color = :black, linewidth = 3)
        lines!(ax, μ3, x3, color = :black, linewidth = 3, linestyle = :dash)
        # Plot the solution up to the current timestep
        lines!(ax, μ[1:n], u[1:n], color = (:teal, 0.75), linewidth = 5)
        # Plot the solution at the current timestep
        scatter!(ax, μ[n], u[n], color = (:teal, 1.0), markersize = 25, strokewidth = 3)

        # Export the phase plane figure 
        save("../fig/slide_04/B-tipping/$n.png", fig)
end

#############
# R-tipping #
#############

# Import the data from csv
u_R = readin("../data/slide_04/R-tipping.csv")
t = u_R[:,1] 
μ = u_R[:,2]
u = u_R[:,3]

# Define the limits for the scalar potential plot
x_inf = -(25.0::Float64)
x_sup = 10.0::Float64
y_inf = -(1500.0::Float64)
y_sup = -y_inf 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = -x*(p^2 - 1.0) - p*x^2 - (1/3)*x^3

# Loop over the parameter values
printstyled("Generating figures for R-tipping using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:brown2)
@showprogress Threads.@threads for n in 1:length(μ)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [800,1200],
                        bg_out = colorant"rgb(147,237,213)",
                        border_color = RGBAf(0,0.236,0.236),
                        box_position = [2:3,1],
                        pad = 30, # Equal amount of "whitespace" in all four directions 
                        limits = ((x_inf, x_sup), (y_inf, y_sup)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[n]) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ[n]), color = (:teal, 1.0), markersize = 45, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_04/R-tipping/$n.png", fig)
end

#############
# N-tipping #
#############

# Import the data from csv
u_N = readin("../data/slide_04/N-tipping.csv")
t = u_N[:,1] 
u = u_N[:,2]

# Define the limits for the scalar potential plot
x_inf = -(1.0::Float64)
x_sup = 3.5::Float64
y_inf = -(0.5::Float64)
y_sup = 2.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Define the fixed parameter value
μ0 = 0.625::Float64

# Loop over the parameter values
printstyled("Generating figures for N-tipping using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:darkgreen)
@showprogress Threads.@threads for n in 1:length(t)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [800,1200],
                        bg_out = colorant"rgb(147,237,213)",
                        border_color = RGBAf(0,0.236,0.236),
                        box_position = [2:3,1],
                        pad = 30, # Equal amount of "whitespace" in all four directions 
                        limits = ((x_inf, x_sup), (y_inf, y_sup)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ0) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ0), color = (:teal, 1.0), markersize = 45, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_04/N-tipping/$n.png", fig)
end
