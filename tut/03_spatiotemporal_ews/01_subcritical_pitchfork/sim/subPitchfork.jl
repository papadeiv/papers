using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets

# Parameters, variables, and derivatives
@parameters t x
@variables u(..)
∂t = Differential(t)
∂x = Differential(x)
∂xx = Differential(x)^2

# Spatio-temporal domain
a = -10.0
b = 10.0
T = 2.0
domains = [t ∈ Interval(0.0, T),
           x ∈ Interval(a, b)]

# Steady-states
u_1 = x -> √(2.0)*sech.(x) 
u_2 = x -> 0.0 .+ 0.0*x
u_3 = x -> 1.0 .+ 0.0*x
u_4 = x -> -1.0 .+ 0.0*x

# ICs
u0_1(t, x) = exp(-10*x^2) 
u0_2(t, x) = (sin(pi*x))^2
u0_3(t, x) = (cos((pi/2)*x))^2
u0_4(t, x) = sech(x)
u0_5(t, x) = 2*sech(x)
u0_6(t, x) = 1 + (1/10)*sech(10*x) 
u0_7(t, x) = 1 - (1/10)*sech(10*x)
u0_8(t, x) = -1 + (1/10)*sech(10*x)
u0_9(t, x) = -1 - (1/10)*sech(10*x)

# IBVP 
pde  = ∂t(u(t, x)) ~ ∂xx(u(t, x)) + (u(t, x))^3 - u(t, x)
bcs = [u(0, x) ~ u0_10(0, x), # IC
       ∂x(u(t, a)) ~ 0.0,     # BC (homogeneous Neumann)
       ∂x(u(t, b)) ~ 0.0]     # BC (homogeneous Neumann)

@named ibvp = PDESystem(pde, bcs, domains, [t, x], [u(t, x)])

# Discretize the PDE problem into a system of ODEs 
dx = 5e-2
order = 2
discretization = MOLFiniteDifference([x => dx], t, approx_order=order)
ode = discretize(ibvp, discretization)

# Solve ODE problem
using OrdinaryDiffEq
sol = solve(ode, Tsit5(), saveat=0.02)

# Plot results and compare with exact solution
x_h = sol[x]
t_h = sol[t]
u_h = sol[u(t, x)]

using CairoMakie
CairoMakie.activate!()
fig1 = Figure(; resolution = (1200, 800), backgroundcolor = :transparent)
ax = Axis(fig1[1, 1],
    #backgroundcolor = :transparent,
    #ylabel = L"λ",
)
lines!(ax, x_h, u_1(x_h), color = :red, linewidth = 4, label=L"u_s^{(1)}(x)=\sqrt{2}\text{sech}(x)")
lines!(ax, x_h, u_2(x_h), color = :blue, linewidth = 4, label=L"u_s^{(2)}(x)=0")
lines!(ax, x_h, u_3(x_h), color = :red, linewidth = 4, label=L"u_s^{(3)}(x)=+1")
lines!(ax, x_h, u_4(x_h), color = :red, linewidth = 4, label=L"u_s^{(4)}(x)=-1")
lines!(ax, x_h, u_h[1, :], color = :black, linewidth = 2, label=L"u_0(x)=u(x,0)")
for i in eachindex(t_h)
    lines!(ax, x_h, u_h[i, :], color = :gray, linewidth = 1, linestyle = :dash)
end
lines!(ax, x_h, u_h[end, :], color = :green, linewidth = 2, label=L"u_T\,(x)=u(x,T\,)")
axislegend(L"u(x,0)=\frac{1}{2}", position = :rb, orientation = :vertical)
save("./test.png", fig1)
