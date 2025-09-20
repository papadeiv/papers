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
$(TYPEDSIGNATURES)

Compute the escape rate of an overdamped particle from a polynomial potential barrier.

The polynomial scalar potential is either explicitly provided by `U::Function` (alongise its second derivative `Uxx::Function`, the location of its local minima `a::Float64` and maxima `b::Float64`) or it is automatically constructed from its monomial `coefficients::Vector{Float64}`, in which case it is assumend to be a cubic.

The escape rate is a large deviation principle dependent on the diffusion level `σ` given to the particle.

## Output
`escape_rate::Tuple`
* `escape_rate.prefactor::Float64`: proportional of the product of the second derivatives evaluated at `a` and `b`; comes from Kramer's formula
* `escape_rate.LDP::Float64`: exponential large deviation principle out of the basin of attraction of `a`

## Example
"""
function estimate_escape(coeffs::Vector, σ)
        # Define a cubic polynomial and its second derivative
        U(x) = coeffs[1]*x + coeffs[2]*(x^2) + coeffs[3]*(x^3)
        Uxx(x) = 2*coeffs[2] + 6*coeffs[3]*x

        # Compute the stable (local minima) and unstable (local maxima) equilibria
        a = +(1/(3*coeffs[3]))*(sqrt((coeffs[2])^2 - 3*coeffs[1]*coeffs[3]) - coeffs[2])
        b = -(1/(3*coeffs[3]))*(sqrt((coeffs[2])^2 - 3*coeffs[1]*coeffs[3]) + coeffs[2])
        
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
