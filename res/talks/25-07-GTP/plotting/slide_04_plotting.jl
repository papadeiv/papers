include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

ews_var = colorant"rgb(136,57,239)"
ews_skw = colorant"rgb(230,69,83)"
CtpFlamingo = colorant"rgb(221, 120, 120)"
CtpYellow = colorant"rgb(229, 200, 144)"
CtpPeach = colorant"rgb(254, 100, 11)"
CtpRed = colorant"rgb(210, 15, 57)"

#############
# B-tipping #
#############

# Import the data from csv
u_B = readin("../data/slide_04/B-tipping.csv")
t = u_B[:,1] 
μ = u_B[:,2]
u = u_B[:,3]

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

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
        fig, ax = mkfig(size = [600,600],
                        bg_out = colorant"rgb(239,241,245)",
                        border_color = colorant"rgb(76,79,105)",
                        pad = 30, # Equal amount of "whitespace" in all four directions 
                        limits = ((x_inf, x_sup), (y_inf, y_sup)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ[n]) for x in domain], color = :black, linewidth = 4.5)
        # Plot the history of the state (ball in the landscape) upto the current timestep
        lines!(ax, u[1:n], [V(u[m],μ[m]) for m in 1:n], color = (colorant"rgb(136,57,239)", 0.50), linewidth = 2)
        # Plot the state (ball in the landscape) at the current timestep
        scatter!(ax, u[n], V(u[n],μ[n]), color = (colorant"rgb(136,57,239)", 1.0), markersize = 45, strokewidth = 3)

        # Export the phase plane figure 
        save("../fig/slide_04/B-tipping/$n.png", fig)
end

#############
# N-tipping #
#############

# Import the data from csv
u_N = readin("../data/slide_04/N-tipping.csv")
t = u_N[:,1] 
u = u_N[:,2]
Nt = length(t)

bins = readin("../data/slide_04/bins.csv")
pdf = readin("../data/slide_04/pdf.csv")
Nb = length(bins)

# Define the limits for the scalar potential plot
x_inf = -(1.0::Float64)
x_sup = 3.5::Float64
y_inf = -(0.1::Float64)
y_sup = 2.25::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define the scalar potential function 
V(x, p) = x*p + x^2 - x^3 + (1/5)*(x^4)

# Define the fixed parameter value
μ0 = 0.625::Float64

# Define the vertical shift to separate the potential from the histogram
shift = 0.85::Float64

# Define a list of plotting indices
idx = collect(StepRange(1, Int64(20), Nt)) 

# Loop over the timesteps 
printstyled("Generating figures for N-tipping using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:darkgreen)
@showprogress Threads.@threads for m in 1:length(idx)
        # Get the current plotting index
        n = idx[m] 

        # Create and customise the scalar potential figure
        fig, ax = mkfig(size = [600,600],
                        bg_out = colorant"rgb(239,241,245)",
                        border_color = colorant"rgb(76,79,105)",
                        pad = 30, # Equal amount of "whitespace" in all four directions 
                        limits = ((x_inf, x_sup), (y_inf, y_sup)),
                        toggle_lab = [false,false],
                        toggle_ticks = [false,false],
                        toggle_ticks_lab = [false,false],
                       )
        # Plot the scalar potential (landscape)
        lines!(ax, domain, [V(x,μ0) + shift for x in domain], color = :black, linewidth = 4.5)
        # Plot the state (ball in the landscape)
        scatter!(ax, u[n], V(u[n],μ0) + shift, color = (ews_skw, 1.0), markersize = 45, strokewidth = 3)
        # Plot the histogram upto the current timestep
        barplot!(ax, bins, pdf[:,n].*(n/Nt), color = pdf[:,n], colormap = [CtpYellow, CtpRed], strokecolor = :black, strokewidth = 1)
        #scatter!(ax, bins, pdf[:,n] .- shift, color = (colorant"rgb(230,69,83)", 0.5), markersize = 25, strokewidth = 1)
        #hist!(ax, u[1:n], bins = bins, offset = -shift, normalization = :pdf, color = :values, strokecolor = :black, strokewidth = 1)
 
        # Export the scalar potential figure 
        save("../fig/slide_04/N-tipping/$m.png", fig)
end
