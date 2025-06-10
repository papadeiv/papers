using ModelingToolkit, DomainSets
using MethodOfLines, OrdinaryDiffEq
using ProgressMeter, Suppressor, LaTeXStrings, CairoMakie

# Parameters, variables, and derivatives
@parameters t x # using ModelingToolkit 
@variables u(..)
Γ = [-1.0, -0.5, 0.02, 0.5]
q = 1.0
α = 2.0
β = 1.0
∂t = Differential(t)
∂x = Differential(x)
∂xx = Differential(x)^2
∂4x = Differential(x)^4

# Spatio-temporal domain
a = -1.0
b = 1.0
T = 1.0
domains = [t ∈ Interval(0.0, T), # using DomainSets
           x ∈ Interval(a, b)]

# Loop over the set of parameters
println("Discretizing and solving the IBVP using ", Threads.nthreads(), " threads\n")
#=Threads.@threads=# for µ in Γ
        # Homogeneous steady-states 
        us0 = 0.0
        us11 = sqrt(Complex((1/(2*β))*(α+sqrt(Complex(α^2 + 4*β*(µ - q^4))))))
        us12 = sqrt(Complex((1/(2*β))*(α-sqrt(Complex(α^2 + 4*β*(µ - q^4))))))
        us21 = -us11
        us22 = -us12

        # Spatially localised perturbations
        u1(x,t) = (1/10)*sech(10*x)
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
        u15(x,t) = sech(10*x + 77) + sech(10*x - 45) - sech(10*x - 100) - sech(10*x + 23) - sech(10*x + 750)

        # Initial condition 
        u0(x,t) = u1(x,t) 

        # Build the IBVP 
        pde = ∂t(u(x,t)) ~ (µ - (q)^2)*u(x,t) - 2*q*∂xx(u(x,t)) - ∂4x(u(x,t)) + α*(u(x,t))^3 - β*(u(x,t))^5
        bcs = [u(x,0) ~ u0(x,0),   # IC
               u(a,t) ~ u(b,t)]    # BC (periodic)

        @named ibvp = PDESystem(pde, bcs, domains, [x,t], [u(x,t)])

        # Discretize the PDE problem into a system of ODEs 
        dx = 1e-1
        order = 2
        discretization = MOLFiniteDifference([x => dx], t, approx_order=order) # using MethodOfLines
        ode = @suppress discretize(ibvp, discretization)

        # Solve ODE problem
        sol = solve(ode, TRBDF2(), saveat=0.1, progress=true) # using OrdinaryDiffEq
        x_h = sol[x]
        t_h = sol[t]
        u_h = sol[u(x,t)]

        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 5)
        fig1 = Figure(; size = (1600, 1066))#, backgroundcolor = :transparent)
        ax = Axis(fig1[1, 1],
                  title = L"u_0(x) = \frac{1}{10}\text{sech}(10x)", 
                #backgroundcolor = :transparent,
                #ylabel = L"λ",
                )
        # Plot the steady-states
        lines!(ax, x_h, us0*ones(length(x_h)), color = :red, linewidth = 2, label=L"u_s^{(0)}(x)=0")
        if imag(us11) == 0
                lines!(ax, x_h, real(us11)*ones(length(x_h)), color = :blue, linewidth = 2, label=L"u_s^{(1,1)}(x)\approx %$(round(real(us11); digits=4))")
                lines!(ax, x_h, real(us21)*ones(length(x_h)), color = :blue, linewidth = 2, label=L"u_s^{(2,1)}(x)=-u_s^{(1,1)}(x)")
        end
        if imag(us12) == 0
                lines!(ax, x_h, real(us12)*ones(length(x_h)), color = :red, linewidth = 2, label=L"u_s^{(1,2)}(x)\approx %$(round(real(us12); digits=4))")
                lines!(ax, x_h, real(us22)*ones(length(x_h)), color = :red, linewidth = 2, label=L"u_s^{(2,2)}(x)=-u_s^{(1,2)}(x)")
        end
        # Plot the transient solution at sampled time instants
        for i in eachindex(t_h)
                lines!(ax, x_h, u_h[:,i], color = :gray, linewidth = 1, linestyle = :dash)
        end
        # Plot the initial perturbation (IC) 
        lines!(ax, x_h, u_h[:,1], color = :black, linewidth = 1, label=L"u_{0}(x)=u(x,t=0)")
        # Plot the asymptotic (t -> ∞) steady-state
        lines!(ax, x_h, u_h[:,end], color = :green, linewidth = 2, label=L"u_T(x)=u(x,t=%$T)")
        # Add the legend
        axislegend(L"\mu=%$µ", position = :rb, orientation = :vertical, labelsize=10.0, titlesize=10.0)
        # Export the figure
        save("../results/pdes/SH35/$µ.png", fig1)
end
