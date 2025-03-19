using Roots, ForwardDiff 

# Find the equilibria and their stability
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
