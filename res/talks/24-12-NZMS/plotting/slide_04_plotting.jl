include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

#############
# B-tipping #
#############

# Import the data from csv
time = readin("../data/slide_04/B-tipping/time.csv")
μ = readin("../data/slide_04/B-tipping/μ.csv")
u = readin("../data/slide_04/B-tipping/u.csv")

# Define the limits for the scalar potential plot
x_inf = -(3.0::Float64)
x_sup = 5.0::Float64
y_inf = -(12.5::Float64)
y_sup = 25.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Loop over the parameter values
printstyled("Generating figures for B-tipping\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1000,1000],
                    pad = 30, # Equal amount of "whitespace" in all four directions 
                    bg_out = "#eeeeeeff",
                    limits = ((x_inf, x_sup), (y_inf, y_sup)),
                    toggle_lab = [false,false],
                    toggle_ticks = [false,false],
                    toggle_ticks_lab = [false,false],
                   )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[n]) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ[n]), color = (:dodgerblue4, 1.0), markersize = 45, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_04/B-tipping/$n.png", fig)
end

#############
# R-tipping #
#############

# Import the data from csv
time = readin("../data/slide_04/R-tipping/time.csv")
μ = readin("../data/slide_04/R-tipping/μ.csv")
u = readin("../data/slide_04/R-tipping/u.csv")

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
printstyled("Generating figures for R-tipping\n"; bold=true, underline=true, color=:brown2)
@showprogress for n in 1:length(μ)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1000,1000],
                    pad = 30, # Equal amount of "whitespace" in all four directions 
                    bg_out = "#eeeeeeff",
                    limits = ((x_inf, x_sup), (y_inf, y_sup)),
                    toggle_lab = [false,false],
                    toggle_ticks = [false,false],
                    toggle_ticks_lab = [false,false],
                   )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[n]) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ[n]), color = (:brown2, 1.0), markersize = 45, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_04/R-tipping/$n.png", fig)
end

#############
# N-tipping #
#############

# Import the data from csv
time = readin("../data/slide_04/N-tipping/time.csv")
u = readin("../data/slide_04/N-tipping/u.csv")

# Define the limits for the scalar potential plot
x_inf = -(1.5::Float64)
x_sup = 3.5::Float64
y_inf = -(1.0::Float64)
y_sup = 3.5::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Define the fixed parameter value
μ0 = 1.0::Float64

# Loop over the parameter values
printstyled("Generating figures for N-tipping\n"; bold=true, underline=true, color=:darkgreen)
@showprogress for n in 1:length(time)
        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1000,1000],
                    pad = 30, # Equal amount of "whitespace" in all four directions 
                    bg_out = "#eeeeeeff",
                    limits = ((x_inf, x_sup), (y_inf, y_sup)),
                    toggle_lab = [false,false],
                    toggle_ticks = [false,false],
                    toggle_ticks_lab = [false,false],
                   )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ0) for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ0), color = (:darkgreen, 1.0), markersize = 45, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_04/N-tipping/$n.png", fig)
end
