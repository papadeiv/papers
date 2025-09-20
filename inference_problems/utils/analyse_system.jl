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
$(TYPEDSIGNATURES)

Numerically computes the zeros of a vector field `f::Function`.

This function can either take a simple function `f(x)` of the state variables or a parametrised vector field `f(x,μ)`, in which case you'll need to provide `μ` as well.

## Keyword arguments
* `domain=[-Inf,Inf]`: interval over which the roots of `f` are sought

## Output
`equilibria::Tuple`
* `equilibria.stable::Vector{Float64}`: stable equilibria of `f`
* `equilibria.unstable::Vector{Float64}`: unstable equilibria of `f`
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
        return (
                stable = stable,
                unstable = unstable
               ) 
end

function get_equilibria(f::Function; domain=[-Inf,Inf])
        # Find the zeros of the function in the specified interval
        equilibria = find_zeros(f, domain[1], domain[2])

        # Define empty arrays to separate stable from unstable equilibria
        stable = Float64[]
        unstable = Float64[]

        # Loop over the equilibria 
        for eq in equilibria
                # Compute the Jacobian at the equilibrium
                stability = ForwardDiff.derivative(f,eq)
                # Determine the stability by checking the sign 
                if stability < 0
                        push!(stable, eq)
                else
                        push!(unstable, eq)
                end
        end

        # Return all the equilibria 
        return (
                stable = stable,
                unstable = unstable
               ) 
end

"""
$(TYPEDSIGNATURES)

Computes the local extrema of a polynomial `V`.

## Keyword Arguments
* `domain::Vector`: the interval over which the extrema are sought (defaults to `[-Inf,Inf]`).

## Output
* `extrema::Vector{Float64}`: stationary points of `V` sorted from smallest to largest 
"""
function get_stationary_points(V::Polynomial; domain=[-Inf,Inf])
        # Differentiate the polynomial
        dV = derivative(V)

        # Find the roots of the first-derivative of the polynomial
        points = roots(dV)

        # Return the points sorted from smallest to largest 
        return sort(points) 
end
