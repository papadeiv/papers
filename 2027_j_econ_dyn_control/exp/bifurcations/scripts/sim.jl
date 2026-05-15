"""
    Simulation script
"""

# Number of iterations of the map
N = convert(Int64, 6e1)

# Fixed parameters
τ = 0.150
ρ = 0.879
α = 0.600
θ = 1.0
a_star = 11.42
e_star = a_star*τ*(1 - ρ) - τ*ρ 
γ = τ + e_star
c_tilde = 1.0
z_tilde = c_tilde/(1 - γ)

# Bifurcation parameters grid
M0 = [[τ, γ, c_tilde],
      [τ, 0.2, c_tilde],
      [τ, 0.1, c_tilde],
      [0.3, γ, c_tilde],
      [0.1, γ, c_tilde],
      [τ, γ, 0.85],
      [τ, γ, 0.5],
      [τ, γ, 0.1]
     ] 

# Initial conditions
g0 = 0.048
A0 = 0.870
e0 = 0.000
L0 = 0.364
U0 = [[g0, A0, e0, L0],             # Red
      [g0, 2*A0, e0, 2*L0],         # Blue
      [g0, 50*A0, e0, 10*L0],       # Peach 
      [g0, 10*A0, e0, L0],          # Mauve
      [g0, 0.85*A0, e0, L0]         # Green
     ]

# Update rule for gt
function g(et, Lt, τ)
        return (et + τ*ρ)*min(θ*Lt, a_star)
end

# Update rule for At
function A(At, et, Lt, τ)
        return (1 + g(et, Lt, τ))*At
end

# Update rule for et
function e(et, Lt, τ)
        return max(0, sqrt(g(et, Lt, τ)*τ*(1 - ρ)) - τ*ρ)
end

# Update rule for Lt
function L(gt, At, et, Lt, μ)
        # Extract the parameters
        τ, γ, c_tilde = μ
        # Compute the per-capita income
        zt = z(gt, At, et, Lt, τ)
        # Define the piece-wise defined function
        if zt >= z_tilde
                return (γ/(τ + e(et, Lt, τ)))*Lt
        elseif zt < z_tilde && zt > c_tilde
                return ((1 - c_tilde/zt)/(τ + e(et, Lt, τ)))*Lt
        else
                return 0
        end
end

# Define the 4-dimensional iterated map
function f(u, μ, n)
        # Extract the bifurcation parameter
        τ, γ, c_tilde = μ

        # Extract the state 
        gt, At, et, Lt = u 

        # Update the state
        g_new = g(et, Lt, τ)
        A_new = A(At, et, Lt, τ)
        e_new = e(et, Lt, τ)
        L_new = L(gt, At, et, Lt, μ)

        # Return the iterated state
        SVector(g_new, A_new, e_new, L_new)
end
