using NonlinearSolve 
using LinearAlgebra
using ForwardDiff 

# Fixed parameters
θ = 0.879 
τ = 0.2 
ρ = 0.5 
α = 0.6 
Ψ = 9.17
γm = 0.85
γw = 0.705

# Bifurcation parameters
mT = 1.0
σI = 1.0
μ = [mT, σI]

# Function to update g
function G(em, ew)
        # Return g_new
        return (0.5*(em + ew) + θ*τ)*Ψ
end

# Implicit function to update em, ew
function E(g_new, Iy)
        T = typeof(g_new)
        e0 = T[0.0, 0.0]

        # Define nonlinear implicit system for em_new, ew_new
        function F(e_new, μ)
                # Extract the state variables
                em_new, ew_new = e_new
                g_new = μ

                # Define the residual surfaces
                F1 = ((em_new + θ*τ)*(em_new + θ*τ + g_new))/((ew_new + θ*τ)*(ew_new + θ*τ + g_new)) - γm/γw
                F2 = ((ew_new + θ*τ)*(ew_new + θ*τ + g_new))/(0.5*(em_new + ew_new) + τ/(2 + Iy)) - (2*γw*g_new)/ρ

                return [F1, F2] 
        end

        return e0
end

function f(u)
        g, em, ew = u
        g_new = G(em, ew)
        hm = (em + θ*τ)/(em + θ*τ + g)
        hw = (ew + θ*τ)/(ew + θ*τ + g)
        Iy = (mT^(1 - α))*(hm/hw)^α - 1
        e_new = E(g_new, Iy)
        return [g_new, e_new[1], e_new[2]]
end

Ueq = [1.612086, 0.0, 0.0]
J = ForwardDiff.jacobian(u -> f(u), Ueq)
if maximum(abs.(eigvals(J))) < 1
        println("\nThe state \n $Ueq \nis a stable fixed point")
else
        println("\nThe state \n $Ueq \nis an unstable fixed point")
end
