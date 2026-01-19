"""
    Simulation script
"""

# Number of iterations of the map
N = convert(Int64, 1e2)

# Global counter for plotting purposes
global idx = 1

# Fixed parameters
θ = 0.879 
τ = 0.2 
ρ = 0.5 
α = 0.6 
Ψ = 9.17
γm = 0.85
γw = 0.705

# Bifurcation parameters
mT = 0.5
σI = 1.5
μ = [mT, σI]

# Function to update g
function G(em, ew)
        # Return g_new
        return (0.5*(em + ew) + θ*τ)*Ψ
end

# Implicit function to update em, ew
function E(g_new, Iy, mT)
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

        # Initialise guesses and solutions of the nonlinear solver
        T = typeof(g_new)
        guesses = [T[0.0, 0.0], T[1.0, 1.0]]
        solutions = Vector{Float64}[]

        # Loop over the guesses
        for e0 in guesses
                # Define the 0-problem of the residual surfaces 
                problem = NonlinearProblem(F, e0, g_new)

                # Solve for and update em_new, ew_new
                solution = solve(problem, TrustRegion())
                push!(solutions, solution.u)
        end

        # Loop over the solutions
        for e_new in solutions
                # Extract the solution for em and ew
                em_new, ew_new = e_new

                # Only return strictly non-negative solutions
                if em_new > 0 && ew_new > 0
                        # Impose the conditions of the optimisation problem
                        if em_new < 1 && ew_new < 1
                                # Return the updated em_new and ew_new
                                return e_new
                        else
                                # Return the updated em_new and ew_new
                                return [0.0, 0.0]
                        end
                end
        end
end

# The update algorithm for em, ew follows from the rules and conditions specified in Lemma 3 (p. ) of Debasis' paper
function update_e(g_new, g, em, ew, mT)
        # Define placeholder variables
        em_new = em 
        ew_new = ew 

        # Compute hm, hw
        hm = (em + θ*τ)/(em + θ*τ + g)
        hw = (ew + θ*τ)/(ew + θ*τ + g)

        # Compute Iy
        Iy = (mT^(1 - α))*(hm/hw)^α - 1

        # Compute the two tresholds
        Em = (2*τ*γm)/(ρ*(2+Iy)) - θ*τ
        gm = ((θ*τ)^2)/Em
        Ew = (2*τ*γw)/(ρ*(2+Iy)) - θ*τ
        gw = ((θ*τ)^2)/Ew

        # Define update algorithm
        if g_new > gm && g_new < gw
                println("case 1\n")
                # Solve the reduced case for em
                em_new = maximum(0.0, 0.5*(γm/ρ - 1)*g_new + sqrt((τ*θ - 0.5*(γm/ρ - 1))^2 + Em*g_new - (θ*τ)^2) - θ*τ)
                ew_new = 0.0 

        elseif g_new > gw && gw > gm
                #=
                # Compute and plot the zero-sets of the implicit system
                zero_sets = zero_set(g_new, Iy)
                plot_zero_set(zero_sets)
                global idx = idx + 1
                =#

                # Solve the implicit nonlinear system to update em and ew
                em_new, ew_new = E(g_new, Iy, mT)

        elseif g_new < gm && gm < gw
                println("case 3\n")
                # Impose the economic constrained of the model (i.e. em, ew strictly non-negative)
                em_new = 0.0
                ew_new = 0.0
        else
                println("case 4: gm = $gm, gw = $gw \n")
        end

        # Return em_new, ew_new
        return em_new, ew_new
end

# Function to update z
function Z(g, em_new, em, ew_new, ew, z, μ)
        # Extract the bifurcation parameters
        mT, σI = μ

        # Compute hm, hw
        hm = (em + θ*τ)/(em + θ*τ + g)
        hw = (ew + θ*τ)/(ew + θ*τ + g)

        # Compute Iy
        Iy = (mT^(1 - α))*(hm/hw)^α - 1

        # Compute the linear (in z) coefficient 
        C = (ρ/(1 + σI/mT + ρ))/(0.5*(em_new + ew_new) + τ/(2 + Iy))

        # Return z_new 
        return C*z 
end

# Function to update the state variables
function update(u, μ)
        # Extract the bifurcation parameters
        mT, σI = μ

        # Extract the state variables at current step n
        g, em, ew, z = u

        # Evaluate the RHS at the current step n
        g_new = G(em, ew)
        em_new, ew_new = update_e(g_new, g, em, ew, mT)
        z_new = Z(g, em_new, em, ew_new, ew, z, μ)

        # Return the updated state at step n+1
        return [g_new, em_new, ew_new, z_new]
end

# Define the (3+1)-dimensional iterated map
function f(u, μ, n)
        # Update the state
        g_new, em_new, ew_new, z_new = update(u, μ)

        # Return the iterated state
        SVector(g_new, em_new, ew_new, z_new)
end
