"""
    Postprocessing script
"""

# Define the domain for F(em_t+1, ew_t+1) = 0
em_range = LinRange(-2.0, 3.5, 3000)
ew_range = LinRange(-2.0, 3.5, 3000)

# Nonlinear system F of surfaces in em_t+1, ew_t+1
function f1(em_new, ew_new, g_new, Iy)
    return ((em_new + θ*τ)*(em_new + θ*τ + g_new))/((ew_new + θ*τ)*(ew_new + θ*τ + g_new)) - γm/γw
end

function f2(em_new, ew_new, g_new, Iy)
    return ((ew_new + θ*τ)*(ew_new + θ*τ + g_new))/(0.5*(em_new + ew_new) + τ/(2 + Iy)) - (2*γw*g_new)/ρ
end

# Address the stiffness of the problem (sharp gradients)
function filter(F1, F2)
        # Create local copies of the surfaces
        Z1 = copy(F1) 
        Z2 = copy(F2)

        # Filter out all the surfaces values above 1 in absolute value
        Z1[abs.(Z1) .> 1] .= NaN
        Z2[abs.(Z2) .> 1] .= NaN

        return Z1, Z2
end

# Construct the zero-sets of F(em_t+1, ew_t+1)
function zero_set(g_new, Iy)
        # Construct the two surfaces z = F(x,y) over the 2-dim domain x -> em_new, y -> ew_new
        z1 = [f1(em_new, ew_new, g_new, Iy) for em_new in em_range, ew_new in ew_range]
        z2 = [f2(em_new, ew_new, g_new, Iy) for em_new in em_range, ew_new in ew_range]

        # Preprocess the surfaces for finding their 0-sets (curves in the (em_t+1,ew_t+1)-plane)
        zero_set_1, zero_set_2 = filter(z1, z2)

        return [zero_set_1, zero_set_2]
end
