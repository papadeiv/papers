"""
    Postprocessing script
"""

# Observable h 
function h(et, gt)
        ht = (et + τ*ρ)/(et + τ*ρ + gt)
        return ht^α 
end

# Observable x 
function x(At, Lt)
        xt = At/Lt
        return xt^(1 - α) 
end

# Observable z 
function z(gt, At, et, Lt)
        # Compute the human capital
        ht = h(et, gt) 
        # Compute the per-capita resources
        xt = x(At, Lt) 
        # Return the per-capita income
        return ht*xt 
end
