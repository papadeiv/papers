"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = (η1 = 3.0, η2 = 1.0, η3 = 0.3)

# Drift terms
f1(x, y, α) = μ.η1 - α - x*(1 + abs(x)) - (1-μ.η3)*y
f2(x, y, α) = α - y*(μ.η3 + abs(x))

# Dynamical system 
function F(z, p)
        # Extract the parameters
        (;η1, η2, η3) = p

        # Extract the coordinates
        x, y = z

        # Define the vector field
        return [η1 - η2 - x*(1 + abs(x)) - (1-η3)*y, η2 - y*(η3 + abs(x))]
end
