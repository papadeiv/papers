# Define Kramer's escape rate
function kramer_escape(U::Function, Uxx::Function, a, b, μ, σ)
        # Compute the pre-factor
        C = (1/(2*pi))*sqrt(abs(Uxx(a, μ))*abs(Uxx(b, μ)))

        # Compute the energy term
        E = 2*((U(b, μ)-U(a, μ))/σ^2)

        # Return the escape rate
        return C*exp(-E)
end
