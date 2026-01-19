# Compute Kramer's escape rate

function kramer_escape(U::Function, Uxx::Function, a, b, σ)
        # Compute the prefactor
        C = (1/(2*pi))*sqrt(Uxx(a)*abs(Uxx(b)))

        # Compute the energy term
        E = 2*((U(b)-U(a))/σ^2)

        # Return the escape rate
        return C*exp(-E)
end
