include("../../../../inc/IO.jl")

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
T = 10.0
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
dx = 10
order = 2
discretization = MOLFiniteDifference([x => dx], t, approx_order=order, grid_align=edge_align)
ode = @suppress discretize(ibvp, discretization)

# Specify the bifurcation location
dt = 0.1
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

# Export the data
writeout(x_h, "../data/figure_05/space.csv")
writeout(t_h, "../data/figure_05/time.csv")
writeout(u_h, "../data/figure_05/solution.csv")
writeout(µ_h, "../data/figure_05/parameter.csv")
writeout(points, "../data/figure_05/random_x.csv")
writeout(timeseries, "../data/figure_05/random_u.csv")

# Execute the postprocessing and plotting scripts
#include("../postprocessing/figure_05_postprocessing.jl")
#include("../plotting/figure_05_plotting.jl")
