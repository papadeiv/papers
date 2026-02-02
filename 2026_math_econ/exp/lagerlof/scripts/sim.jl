"""
    Simulation script
"""

# Number of iterations of the map
N = convert(Int64, 3e4)

# Fixed parameters
τ = 0.50
ρ = 0.879
α = 0.6
γ = 0.225
θ = 1.0
a_star = 11.42
c_tilde = 1.0
z_tilde = c_tilde/(1 - γ)
μ = 1.0

# Initial condition
g0 = 0.048
A0 = 0.870
e0 = 0.000
L0 = 0.364
u0 = [g0, A0, e0, L0]

# Update rule for gt
function g(et, Lt)
        return (et + τ*ρ)*min(θ*Lt, a_star)
end

# Update rule for At
function A(At, et, Lt)
        return (1 + g(et, Lt))*At
end

# Update rule for et
function e(et, Lt)
        return max(0, sqrt(g(et, Lt)*τ*(1 - ρ)) - τ*ρ)
end

# Update rule for Lt
function L(gt, At, et, Lt)
        #println("    sqrt(g(et,Lt)τ(1-ρ)) - τρ = $(sqrt(g(et, Lt)*τ*(1 - ρ)) - τ*ρ)")
        # Compute the per-capita income
        zt = z(gt, At, et, Lt)
        #println("    zt = $zt")
        # Define the piece-wise defined function
        if zt >= z_tilde
                #println("    zt > z*")
                return (γ/(τ + e(et, Lt)))*Lt
        elseif zt < z_tilde && zt > c_tilde
                #println("    c* < zt < z*")
                return ((1 - c_tilde/zt)/(τ + e(et, Lt)))*Lt
        else
                #println("    zt < c*")
                return 0
        end
end

# Define the 4-dimensional iterated map
function f(u, μ, n)
        # Extract the bifurcation parameters
        # p, q = μ
        
        # Extract the state 
        gt, At, et, Lt = u 
        #println("gt = $gt, At = $At, et = $et, Lt = $Lt")
        #println("    min(θLt, a*) = $(min(θ*Lt, a_star))")

        # Update the state
        g_new = g(et, Lt) 
        A_new = A(At, et, Lt)
        e_new = e(et, Lt)
        L_new = L(gt, At, et, Lt)
        #println("gt+1 = $g_new, At+1 = $A_new, et+1 = $e_new, Lt+1 = $L_new \n")

        # Return the iterated state
        SVector(g_new, A_new, e_new, L_new)
end
