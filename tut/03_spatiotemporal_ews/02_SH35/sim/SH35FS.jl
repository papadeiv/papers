###############################################################
#                                                             #
#    Cubic-Quintic Swift-Hohenberg with Fast-Slow dynamics    #
#                                                             #
###############################################################

# Initialise a logging terminal session
using Logging:global_logger
using TerminalLoggers:TerminalLogger
global_logger(TerminalLogger())

# Variables and unkowns of the PDE
using ModelingToolkit
@parameters t x 
@variables u(..) µ(..)

# IC of the parameter
µs0 = -0.05

# Piece-wise defined function for the (homogeneous in space) ramping of μ(t)
ε = 0.001
f(t) = t <= 1000 ? ε : 0.0 # Ternary operator
@register_symbolic f(t)

# Parameters of the PDE (non-bifurcating ones)
q = 1.0
α = 2.0
β = 1.0

# Differential operators
∂t = Differential(t)
∂x = Differential(x)
∂xx = Differential(x)^2
∂4x = Differential(x)^4

# Spatio-temporal domain
using DomainSets
a = -32.0*π
b = 32.0*π
T = 400.0
domains = [t ∈ Interval(0.0, T),
           x ∈ Interval(a, b)]

# Homogeneous steady-states 
us0 = 0.0
us11 = sqrt(Complex((1/(2*β))*(α+sqrt(Complex(α^2 + 4*β*(µs0 - q^4))))))
us12 = sqrt(Complex((1/(2*β))*(α-sqrt(Complex(α^2 + 4*β*(µs0 - q^4))))))
us21 = -us11
us22 = -us12

# Spatially localised perturbations
u1(x,t) = sech(x)
u2(x,t) = (1/10)*exp(-10*x^2)
# Spatially periodic perturbations (Lc = 2π)
u3(x,t) = (1/10)*sin(x)
u4(x,t) = (1/10)*sin(x) + us11
u4_2(x,t) = sin(x) + us11
u5(x,t) = sinh(cos(x))
u6(x,t) = exp(cos(x))
# Spatially periodic perturbations (L = π < Lc)
u7(x,t) = sech(cos(x))
u8(x,t) = sech(cos(x)) + us11
u9(x,t) = sech(x*sin(x))
# Flat perturbations
u10(x,t) = us0 + 1e-1 
u11(x,t) = us0 - 1e-1
u12(x,t) = us11 + 1e-1
u13(x,t) = us11 - 1e-1
u14(x,t) = us12 - 1e-1
# Multiple compact pulses
u15(x,t) = sech(10*x + 570) + sech(10*x - 450) - sech(10*x - 100) - sech(10*x + 230) - sech(10*x + 750)

# Initial conditions
u0(x,t) = u15(x,t) 
µ0(x,t) = µs0

# Build the IBVP 
pde = [∂t(u(x,t)) ~ (µ(x,t) - (q)^2)*u(x,t) - 2*q*∂xx(u(x,t)) - ∂4x(u(x,t)) + α*(u(x,t))^3 - β*(u(x,t))^5,
       ∂t(µ(x,t)) ~ f(t)]
bcs = [# Solution
       u(x,0) ~ u0(x,0),   # IC
       u(a,t) ~ u(b,t),    # BC (periodic)
       # Parameter
       µ(x,0) ~ µ0(x,0),   # IC
       µ(a,t) ~ µ(b,t)     # BC (periodic)
      ]
@named ibvp = PDESystem(pde, bcs, domains, [x,t], [u(x,t), µ(x,t)])

# Discretize the PDE problem into a system of ODEs 
using MethodOfLines, Suppressor
println("Discretizing and solving the IBVP using ", Threads.nthreads(), " threads")
dx = 1000 
order = 2
discretization = MOLFiniteDifference([x => dx], t, approx_order=order, grid_align=edge_align)
ode = @suppress discretize(ibvp, discretization)

# Specify the bifurcation location
dt = 0.001
µ_crit = 0.0 
T_crit = (µ_crit-µs0)/ε
dt_crit = trunc(Int, T_crit/dt) 

# Solve the discretised ODE problem
using OrdinaryDiffEq
sol = solve(ode, TRBDF2(), saveat=dt, progress=true)

# Extract the data
x_h = sol[x]
t_h = sol[t]
u_h = sol[u(x,t)]
µ_h = sol[µ(x,t)]

# Extract the timeseries of the solution value at random points in the domain
using StatsBase
points = sample(1:size(x_h,1), 9; replace=false)
timeseries = [u_h[x,:] for x in points] 

# Plot the timeseries
using CairoMakie, Makie.Colors, LaTeXStrings
CairoMakie.activate!(; px_per_unit = 5)
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
idx = CartesianIndices((1:3, 1:3))
for i ∈ eachindex(timeseries)
    pos = [idx[i][1], idx[i][2]] 
    local ax = Axis(fig[idx[i][1], idx[i][2]],
              backgroundcolor = :transparent,
              #title = L"u(x=%$(x_h[points[i]]),t)",
              #xlabel = L"t"
              xgridvisible = false,
              xticksize = 10,
              xticklabelsize = 24,
              ygridvisible = false,
              yticksize = 10,
              yticklabelsize = 24,
    )
    lines!(ax, t_h, timeseries[i], color = i, colormap = :tab10, colorrange = (1,9), linewidth = 2)
    # Plot a vertical line signaling the bifurcation of the homogeneous steady state
    lines!(ax, T_crit.*ones(2), [timeseries[i][argmin(timeseries[i])], timeseries[i][argmax(timeseries[i])]], color = :black, linestyle = :dash)
    # Export the figure
    if i == size(timeseries, 1)
            save("../results/pdes/SH35FS/fig3.5.2.png", fig)
    end
end

# Compute the norms of the solution at each timestep
using LinearAlgebra
u_L1 = (1/(abs(b-a))).*[norm(u_h[:,t], 1) for t ∈ eachindex(t_h)]
#u_L2 = (1/(abs(b-a))).*[norm(u_h[:,t], 2) for t ∈ eachindex(t_h)]
u_L2 = [norm(u_h[:,t], 2) for t ∈ eachindex(t_h)]
u_L∞ = [norm(u_h[:,t], Inf) for t ∈ eachindex(t_h)]

# Export the L2 norm to csv
using CSV, Tables
mat = Matrix{Float64}(undef, length(u_L2), 2)
mat[:,1] = µ_h[1,:]
mat[:,2] = u_L2
CSV.write("./pde2path/sh35/solution.csv", Tables.table(mat), delim=',', writeheader=false)

# Dynamic Mode Decomposition
using Printf
@printf("EWS 1: leading eigenvalue of the DMD approximation of the Koopman operator\n")

window_size = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50]
window_stride = 1

using ProgressMeter
CairoMakie.activate!(; px_per_unit = 2)
@showprogress for width in window_size
        # Eigenvalues of the DMD (approximate spectrum of the Koopman operator)
        λ = zeros(Float64, width, (size(t_h, 1) - width))
        λmax = zeros(Float64, (size(t_h, 1) - width))
        # Eigenvalues of the SVD (singular values)
        σ = zeros(Float64, width, (size(t_h, 1) - width))
        σmax = zeros(Float64, (size(t_h, 1) - width))

        Threads.@threads for k=1:(size(t_h, 1) - width)
                # Assemble the first snapshot matrix
                X = hcat([u_h[:, t] for t ∈ k:k+width-1]...) #display(X)
                # Perform the SVD
                F = svd(X)
                # Extract the singular values
                σ[:,k] = real(F.S)
                σmax[k] = σ[1,k]
                # Assemble the second snapshot matrix
                Y = hcat([u_h[:, t] for t ∈ (k+1):k+width]...)
                # Perform the DMD 
                S = (F.U)'*Y*(F.Vt)'*inv(diagm(F.S)) 
                Λ = eigvals(S)
                Q = eigvecs(S)
                # Extract the eigenvalues 
                λ[:,k] = sort(real(Λ), rev=true)
                λmax[k] = λ[1,k]
                # Plot both decays at current timestep
                #=
                local fig = Figure(; size = (1200, 800))
                local ax1 = Axis(fig[1, 1], yticklabelcolor = :blue, yscale=log10)
                local ax2 = Axis(fig[1, 1], yticklabelcolor = :red, yaxisposition = :right)
                hidespines!(ax2)
                hidexdecorations!(ax2)
                local sings = scatterlines!(ax1, 1:width, σ[:,k], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 14)
                local eigns = scatterlines!(ax2, 1:width, λ[:,k], color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 14)
                # Add the legend
                axislegend(ax2, [sings, eigns], [L"\text{singular values of } \mathcal{X}", L"\text{eigenvalues of } \mathcal{K}"], "t=$(round(t_h[k]; digits = 4))", position = :rt)
                # Export the image
                save("../results/pdes/SH35FS/DMD/decay/$width/$k.png", fig)
                =#
        end

        # Plot the early-warning signal from the DMD
        cap = 1.002
        shoe = 0.9995
        ƛ = λmax[shoe .< λmax .< cap]
        local fig = Figure(; size = (1200, 400), backgroundcolor = :transparent)
        local ax = Axis(fig[1, 1],
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((0, T), nothing),#(shoe,cap)),
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabel = L"\text{time}",
                xlabelvisible = true,
                xlabelsize = 32,
                xlabelcolor = :black,
                xlabelpadding = -20.0,
                xticks = [0, T],#t_h[size(ƛ,1)]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 24,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.0f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabel = L"\lambda_{\text{max}}",
                ylabelvisible = true,
                ylabelsize = 32,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [0,2.5],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 24,
                ytickformat = "{:.3f}",
                yscale = identity,
                yaxisposition = :left,
                yticklabelcolor = :black,
        )
        lines!(ax, LinRange(t_h[width], t_h[dt_crit], size(λmax[100:dt_crit], 1)), λmax[100:dt_crit]) # label = "Accuracy = $(round((size(ƛ,1)/size(λmax,1))*100; digits = 1))%")
        # Plot a vertical line signaling the bifurcation of the homogeneous steady state
        lines!(ax, T_crit.*ones(2), [λmax[argmin(λmax[100:dt_crit])], λmax[argmax(λmax[100:dt_crit])]], color = :black, linestyle = :dash)
        # Export the figure
        save("../results/pdes/SH35FS/DMD/eigenvalues/$width.png", fig)

        # Plot the early-warning signal from the SVD
        local fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
        # Plot the full (temporal) range of the leading eigenvalue
        local ax1 = Axis(fig[1, 1:2],
                xlabel = L"\text{time}",
                ylabel = L"\sigma_{\text{max}}",
                limits = (0, T, nothing, nothing),
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabelvisible = true,
                xlabelsize = 24,
                xlabelcolor = :black,
                xlabelpadding = -20.0,
                xticks = [0, T],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 24,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.0f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabelvisible = true,
                ylabelsize = 28,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                #yticks = [σmax[argmin(σmax)], σmax[argmax(σmax)]],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 24,
                ytickformat = "{:.0f}",
                yscale = identity,
                yaxisposition = :left,
                yticklabelcolor = :black,
        )
        lines!(ax1, LinRange(t_h[width], t_h[end], size(σmax,1)), σmax, label = "Leading singular values")
        # Plot a vertical line signaling the bifurcation of the homogeneous steady state
        lines!(ax1, T_crit.*ones(2), [σmax[argmin(σmax)], σmax[argmax(σmax)]], color = :black, linestyle = :dash)
        # Extract a detailed temporal range for the EWS
        #=
        inf_t = 400
        sup_t = 1300
        σres = σmax[inf_t:sup_t]
        # Plot a red box in the previous figure showing the detailed temporal evolution in the bottom subplot that follows
        using Makie.GeometryBasics
        poly!(ax1, Point2f[(t_h[inf_t], -10), (t_h[sup_t], -10), (t_h[sup_t], 10), (t_h[inf_t], 10)], color = (:white, 0.0), strokecolor = :black, strokewidth = 1)
        # Plot a detailed range of the temporal evolution of the leading eigenvalue
        local ax2 = Axis(fig[2, 1:2],
                xlabel = L"\text{time}",
                ylabel = L"\sigma_{\text{max}}",
                limits = (t_h[inf_t], t_h[sup_t], nothing, nothing),
                # Background
                backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                # Title
                #title = (L"\mu=%$(round(Y; digits=3))"),
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # x-axis
                xlabelvisible = true,
                xlabelsize = 24,
                xlabelcolor = :black,
                xlabelpadding = -20.0,
                xticks = [t_h[inf_t],t_h[sup_t]],
                xticksvisible = true,
                xticksize = 10,
                xticklabelsvisible = true,
                xticklabelsize = 24,
                xticklabelalign = (:center,:top),
                xtickformat = "{:.0f}",
                xscale = identity, #log10,
                xaxisposition = :bottom,
                # y-axis
                ylabelvisible = true,
                ylabelsize = 28,
                ylabelcolor = :black,
                ylabelpadding = -30.0,
                yticks = [σres[argmin(σres)], σres[argmax(σres)]],
                yticksvisible = true,
                yticksize = 10,
                yticklabelsvisible = true,
                yticklabelsize = 24,
                ytickformat = "{:.2f}",
                yscale = identity,
                yaxisposition = :left,
                yticklabelcolor = :black,
        )
        lines!(ax2, t_h[inf_t:sup_t], σres)
        # Plot a vertical line signaling the bifurcation of the homogeneous steady state
        lines!(ax2, 9.0.*ones(2), [σres[argmin(σres)], σres[argmax(σres)]], color = :black, linestyle = :dash)
        # Plot a vertical line signaling the visual detection of the pattern formation
        #lines!(ax, 16.0.*ones(2), [shoe, cap*2]#=[σmax[argmin(σmax)], σmax[argmax(σmax)]]=#, color = :black, linestyle = :dash)
        # Plot a vertical line signaling the end of the parameter's ramping
        # Plot a vertical line signaling the reaching of the steady-state solution 
        #lines!(ax1, 30.0.*ones(2), [σres[argmin(σres)], σres[argmax(σres)]], color = :red)
        # Add the legend
        #axislegend("Window size = $width snapshots", position = :lt)
        =#
        # Export the figure
        save("../results/pdes/SH35FS/DMD/singularvalues/$width.png", fig)
end

# Spatial variance
@printf("EWS 2: spatial variance of random sampling nodes\n")
sampling_nodes = [3, 9, 27, 81, 243, 729] 
@showprogress for samples ∈ sampling_nodes
        # Extract random sampling nodes from the domain
        nodes = sample(1:size(x_h,1), samples; replace=false)
        # Extract the timeseries of the solution value at each sampling node
        local timeseries = hcat([u_h[node,:] for node ∈ nodes]...)
        # Define spatial mean and variance 
        local M = zeros(Float64, size(t_h, 1))
        local V = zeros(Float64, size(t_h, 1))
        Threads.@threads for t ∈ eachindex(t_h)
                # Compute the spatial mean among the sampling nodes at timestep t
                M[t] = mean(timeseries[t,:]) 
                # Compute the spatial variance of the sampling nodes at timestep t
                V[t] = var(timeseries[t,:]; mean=M[t]) 
        end

        # Create the plot on a twin axis 
        local fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
        local ax1 = Axis(fig[1, 1], 
              # Background
              backgroundcolor = :transparent,
              xgridvisible = false,
              ygridvisible = false,
              limits = ((0, T), nothing),
              # Title
              #title = (L"\mu=%$(round(Y; digits=3))"),
              titlevisible = false,
              titlesize = 22,
              titlealign = :center,
              titlegap = 4.0,
              # x-axis
              xlabel = L"\text{time}",
              xlabelvisible = true,
              xlabelsize = 24,
              xlabelcolor = :black,
              xlabelpadding = -20.0,
              xticks = [0, T],
              xticksvisible = true,
              xticksize = 10,
              xticklabelsvisible = true,
              xticklabelsize = 24,
              xticklabelalign = (:center,:top),
              xtickformat = "{:.0f}",
              xscale = identity, #log10,
              xaxisposition = :bottom,
              # y-axis
              ylabel = L"\text{spatial mean}",
              ylabelvisible = true,
              ylabelsize = 36,
              ylabelcolor = :black,
              ylabelpadding = -30.0,
              yticks = [-0.01,0],
              yticksvisible = true,
              yticksize = 10,
              yticklabelsvisible = true,
              yticklabelsize = 24,
              ytickformat = "{:.1f}",
              yscale = identity,
              yaxisposition = :left,
              yticklabelcolor = :blue,
        )
        local ax2 = Axis(fig[1, 1], 
              # Background
              backgroundcolor = :transparent,
              xgridvisible = false,
              ygridvisible = false,
              limits = ((0, T), nothing),
              # Title
              #title = (L"\mu=%$(round(Y; digits=3))"),
              titlevisible = false,
              titlesize = 22,
              titlealign = :center,
              titlegap = 4.0,
              # x-axis
              xlabel = L"\text{time}",
              xlabelvisible = false,
              xlabelsize = 24,
              xlabelcolor = :black,
              xlabelpadding = -20.0,
              xticks = [0.5, 0.68, 1.2],
              xticksvisible = false,
              xticksize = 10,
              xticklabelsvisible = false,
              xticklabelsize = 24,
              xticklabelalign = (:right,:top),
              xtickformat = "{:.2f}",
              xscale = identity, #log10,
              xaxisposition = :bottom,
              # y-axis
              ylabel = L"\text{spatial variance}",
              ylabelvisible = true,
              ylabelsize = 24,
              ylabelcolor = :black,
              ylabelpadding = -30.0,
              yticks = [0,1.5],
              yticksvisible = true,
              yticksize = 10,
              yticklabelsvisible = true,
              yticklabelsize = 24,
              ytickformat = "{:.1f}",
              yscale = identity,
              yaxisposition = :right,
              yticklabelcolor = :red,
        )
        hidespines!(ax2)
        hidexdecorations!(ax2)
        local meanline = lines!(ax1, t_h, M, color = :blue, label = L"\mathbb{E}(f(x_k))")
        local varline = lines!(ax2, t_h, V, color = :red, label = L"\text{Var}(f(x_k))")
        # Plot a vertical line signaling the bifurcation of the homogeneous steady state
        lines!(ax1, T_crit.*ones(2), [M[argmin(M)], M[argmax(M)]], color = :black, linestyle = :dash)
        # Add the legend
        #axislegend(ax2, [meanline, varline], [L"\mathbb{E}(f(x_k))", L"\text{Var}(f(x_k))"], "Spatial indicators", position = :lt)
        # Export the image
        save("../results/pdes/SH35FS/SV/$samples.png", fig)
end

# Compute and plot the spatial indicators above but for all the nodes in the discretised domain
M = [mean(u_h[:,t]) for t ∈ eachindex(t_h)]
V = [var(u_h[:,t]; mean=M[t]) for t ∈ eachindex(t_h)]
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
ax1 = Axis(fig[1, 1], 
              # Background
              backgroundcolor = :transparent,
              xgridvisible = false,
              ygridvisible = false,
              limits = ((0, T), nothing),
              # Title
              #title = (L"\mu=%$(round(Y; digits=3))"),
              titlevisible = false,
              titlesize = 22,
              titlealign = :center,
              titlegap = 4.0,
              # x-axis
              xlabel = L"\text{time}",
              xlabelvisible = true,
              xlabelsize = 40,
              xlabelcolor = :black,
              xlabelpadding = -20.0,
              xticks = [0, T],
              xticksvisible = true,
              xticksize = 10,
              xticklabelsvisible = true,
              xticklabelsize = 28,
              xticklabelalign = (:center,:top),
              xtickformat = "{:.0f}",
              xscale = identity, #log10,
              xaxisposition = :bottom,
              # y-axis
              ylabel = L"\text{spatial mean}",
              ylabelvisible = true,
              ylabelsize = 40,
              ylabelcolor = :blue,
              ylabelpadding = -30.0,
              yticks = [M[argmin(M)],0],
              yticksvisible = true,
              yticksize = 10,
              yticklabelsvisible = true,
              yticklabelsize = 28,
              ytickformat = "{:.2f}",
              yscale = identity,
              yaxisposition = :left,
              yticklabelcolor = :blue,
)
ax2 = Axis(fig[1, 1], 
              # Background
              backgroundcolor = :transparent,
              xgridvisible = false,
              ygridvisible = false,
              limits = ((0, T), nothing),
              # Title
              #title = (L"\mu=%$(round(Y; digits=3))"),
              titlevisible = false,
              titlesize = 22,
              titlealign = :center,
              titlegap = 4.0,
              # x-axis
              xlabel = L"\text{time}",
              xlabelvisible = false,
              xlabelsize = 24,
              xlabelcolor = :black,
              xlabelpadding = -20.0,
              xticks = [0.5, 0.68, 1.2],
              xticksvisible = false,
              xticksize = 10,
              xticklabelsvisible = false,
              xticklabelsize = 28,
              xticklabelalign = (:right,:top),
              xtickformat = "{:.2f}",
              xscale = identity, #log10,
              xaxisposition = :bottom,
              # y-axis
              ylabel = L"\text{spatial variance}",
              ylabelvisible = true,
              ylabelsize = 40,
              ylabelcolor = :red,
              ylabelpadding = -30.0,
              yticks = [0,V[argmax(V)]],
              yticksvisible = true,
              yticksize = 10,
              yticklabelsvisible = true,
              yticklabelsize = 28,
              ytickformat = "{:.1f}",
              yscale = identity,
              yaxisposition = :right,
              yticklabelcolor = :red,
)
hidespines!(ax2)
hidexdecorations!(ax2)
meanline = lines!(ax1, t_h, M, color = :blue, label = L"\mathbb{E}(f(x_k))")
varline = lines!(ax2, t_h, V, color = :red, label = L"\text{Var}(f(x_k))")
# Plot a vertical line signaling the bifurcation of the homogeneous steady state
lines!(ax1, T_crit.*ones(2), [M[argmin(M)], M[argmax(M)]], color = :black, linestyle = :dash)
# Add the legend
#axislegend(ax2, [meanline, varline], [L"\mathbb{E}(f(x_k))", L"\text{Var}(f(x_k))"], "Spatial indicators", position = :lt)
# Export the image
save("../results/pdes/SH35FS/fig3.5.1.png", fig)

# Plot of the time-evolution of the solution
@printf("Exporting the plots of the %d snapshots", size(t_h, 1))
CairoMakie.activate!(; px_per_unit = 2)
Threads.@threads for i ∈ eachindex(t_h)
        fig = Figure(; size = (1200, 800), backgroundcolor = "#fdfff2ff")
        ax = Axis(fig[1, 1],
                  backgroundcolor = "#fdfff2ff",
                  #title = latexstring("u_0(x) = \\text{Train of pulses}, \\quad t=", @sprintf("%.3f", t_h[i]), ", \\quad µ =", @sprintf("%.3f", µ_h[1,i])),
                  xticksvisible = false,
                  #xticklabelsvisible = false,
                  xgridvisible = false,
                  yticksize = 10,
                  yticklabelsize = 24,
                  xticksize = 10,
                  xticklabelsize = 24,
                  yticksvisible = false,
                  #yticklabelsvisible = false,
                  ygridvisible = false,
                  limits = (a, b, -2.0, 2.0),
                  xlabel = L"x",
                  xlabelsize = 32,
                  xlabelcolor = :black,
                  xlabelpadding = -20.0,
                  xticks = [-100,-50,50,100],
                  ylabel = L"u",
                  ylabelsize = 32,
                  ylabelcolor = :black,
                  ylabelpadding = -20.0,
                  yticks = [-2,-1,1,2],
        )
        # Plot the steady-states
        #lines!(ax, x_h, us0*ones(length(x_h)), color = :green, linewidth = 2, label=L"u_s^{(0)}(x)=0")
        if imag(us11) == 0
                lines!(ax, x_h, real(us11)*ones(length(x_h)), color = :blue, linewidth = 2, label=L"u_s^{(1,1)}(x)\approx %$(round(real(us11); digits=4))")
                lines!(ax, x_h, real(us21)*ones(length(x_h)), color = :blue, linewidth = 2, label=L"u_s^{(2,1)}(x)=-u_s^{(1,1)}(x)")
        end
        if imag(us12) == 0
                lines!(ax, x_h, real(us12)*ones(length(x_h)), color = :red, linewidth = 2, label=L"u_s^{(1,2)}(x)\approx %$(round(real(us12); digits=4))")
                lines!(ax, x_h, real(us22)*ones(length(x_h)), color = :red, linewidth = 2, label=L"u_s^{(2,2)}(x)=-u_s^{(1,2)}(x)")
        end
        # Plot the transient solution at current time instant
        lines!(ax, x_h, u_h[:,i], color = :black, linewidth = 2)
        #text!(0, 1.75, text = L"\text{(c)}", align = [:center, :center], color = :black, fontsize = 24)
        # Plot the initial perturbation (IC) 
        #lines!(ax, x_h, u_h[:,1], color = :orange, linewidth = 1, label=L"u_{0}(x)=u(x,t=0)")
        # Plot the random points extracted earlier
        for j ∈ eachindex(timeseries)
                #lines!(ax, x_h[points[j]].*ones(2), [min(0.0, timeseries[j][i]), max(0.0, timeseries[j][i])], color = j, colormap = :tab10, colorrange = (1,9), linestyle = :dash)
                #scatter!(ax, x_h[points[j]], timeseries[j][i], color = j, colormap = :tab10, colorrange = (1,9), strokecolor = :black, strokewidth = 1.5, markersize = 15)
        end
        # Export the figure
        save("../results/pdes/SH35FS/flipbook/$i.png", fig)
end
