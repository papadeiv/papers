using ApproxFun, LinearAlgebra

# Spatio-temporal domain
a = -10
b = 10
dx = a..b
dt = 0..2
Ω = dx × dt

# Differential operators
∂t = Derivative(Ω,[0,1])
∂x = Derivative(Ω,[1,0])

# ICs
u1 = Fun(x->sech(10*x), dx)
u2 = Fun(x->exp(-10*x^2), dx)
u3 = Fun(x->sinpi(x), dx)
u4 = Fun(x->cospi(x/2), dx)

# Linear operator
L = [∂t - ∂x^2 + 1; I⊗ldirichlet(dt); bvp(dx)⊗I];

# RHS with IC and BCs
RHS = [0.0; u2; 0.0; 0.0]

# Forward-time evolution
u = \(L, RHS, tolerance=1e-2)

# Plots
import Plots
xplot = a:0.02:b
p = Plots.plot(xplot, u.(xplot, 0), label="t=0", legend=true, linewidth=2)
for t in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
	Plots.plot!(xplot, u.(xplot, t), label="t=$t")
end
p
