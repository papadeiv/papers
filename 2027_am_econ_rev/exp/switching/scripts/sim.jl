"""
    Simulation script
"""

# Number of iterations of the map
N = convert(Int64, 3e1)

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
μ = 1.0

# Initial conditions
g0 = 0.048
A0 = 0.870
e0 = 0.000
L0 = 0.364
U0 = [#[g0, A0, e0, L0],       # Red  (Lagerlof's)
      [g0, 5, e0, 0.5],        # Blue  (white region)
      [g0, 5, e0, 3.5],        # Green (blue region)
      [g0, 5, e0, 1.5],        # Peach (pink region) 
      [g0, 5, e0, 3.0],        # Mauve (white region)
     ]

# Update rule for gt
function g(et, Lt)
        println("* Update for g: ")
        if Lt < a_star
                println("    Lt < a_star")
        else
                println("    Lt ≥ a_star")
        end

        return (et + τ*ρ)*min(θ*Lt, a_star)
end

# Update rule for At
function A(At, et, Lt)
        return (1 + g(et, Lt))*At
end

# Update rule for et
function e(et, Lt)
        println("* Update for e: ")
        if sqrt(g(et, Lt)*τ*(1 - ρ)) < τ*ρ
                println("    sqrt() - τρ < 0")
        else
                println("    sqrt() - τρ ≥ 0")
        end
        return max(0, sqrt(g(et, Lt)*τ*(1 - ρ)) - τ*ρ)
end

# Update rule for Lt
function L(gt, At, et, Lt)
        println("* Update for L: ")
        # Compute the per-capita income
        zt = z(gt, At, et, Lt)
        # Define the piece-wise defined function
        if zt >= z_tilde
                println("    z = $zt >= z_tilde = $z_tilde")
                return (γ/(τ + e(et, Lt)))*Lt
        elseif zt < z_tilde && zt > c_tilde
                println("    c_tilde = $c_tilde <= z = $zt < z_tilde = $z_tilde")
                return ((1 - c_tilde/zt)/(τ + e(et, Lt)))*Lt
        else
                println("    z = $zt < c_tilde = $c_tilde")
                return 0
        end
end

# Define the 4-dimensional iterated map
function f(u, μ, n)
        # Extract the state 
        gt, At, et, Lt = u 
        println("At = $At, Lt = $Lt, gt = $gt, et = $et")

        # Update the state
        g_new = g(et, Lt) 
        A_new = A(At, et, Lt)
        e_new = e(et, Lt)
        L_new = L(gt, At, et, Lt)
        println("At = $A_new, Lt = $L_new, gt = $g_new, et = $e_new")

        # Return the iterated state
        SVector(g_new, A_new, e_new, L_new)
end
