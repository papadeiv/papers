"""
Function utilities for the application of Freidlin-Wentzell theory to bounded sample paths. 

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#----------------------------------#
#                                  # 
#   large_deviation_principle.jl   #                     
#                                  #
#----------------------------------#

"""
    estimate_escape(U, a, b, σ; kwargs...)

Compute the LDP to escape the basin of attraction of a potential `U` for a particle starting at the minimum `a` and crossing a maximum `b` subject to diffusion level `σ`.
"""
function estimate_escape(U::Function, Uxx::Function, a, b, σ)
        # Compute the prefactor
        C = (1/(2*pi))*sqrt(Uxx(a)*abs(Uxx(b)))

        # Compute the energy term
        E = 2*((U(b)-U(a))/σ^2)

        # Return the escape rate
        return (
                prefactor = C,
                LDP = exp(-E)
               )
end
