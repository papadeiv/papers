"""
Dynamics utility functions to compute the typical analytical properties of dynamical systems (equilibria, spectrum etc...).

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#-----------------------#
#                       # 
#   analyse_system.jl   #                     
#                       #
#-----------------------#

"""
    get_equilibria(f, μ; domain)

Given a vector field `f` and its bifurcation parameter `μ` find all its equilibria and their stability.
Detailed description.
"""
function get_equilibria(f::Function, μ::Float64; domain=[-Inf,Inf])
        # Reduce the parameter dependent dynamics to a 1D scalar function
        F(x) = f(x, μ)

        # Find the zeros of the function in the specified interval
        equilibria = find_zeros(F, domain[1], domain[2])

        # Define empty arrays to separate stable from unstable equilibria
        stable = Float64[]
        unstable = Float64[]

        # Loop over the equilibria 
        for eq in equilibria
                # Compute the Jacobian at the equilibrium
                stability = ForwardDiff.derivative(F,eq)
                # Determine the stability by checking the sign 
                if stability < 0
                        push!(stable, eq)
                else
                        push!(unstable, eq)
                end
        end

        # Return all the equilibria 
        return [stable, unstable] 
end
