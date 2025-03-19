# Approximate the potential function from the data using polynomial interpolation
function approximate_potential(data_x, data_f; degree=3::Int64)
        interpolant = Polynomials.fit(data_x, data_f, degree)
        return interpolant
end
