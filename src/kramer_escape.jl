# Define Kramer's escape rate
function kramer_escape(U::Function, Uxx::Function, a, b, σ)
        # Compute the pre-factor
        C = (1/(2*pi))*sqrt(abs(Uxx(a))*abs(Uxx(b)))

        # Compute the energy term
        E = 2*((U(b)-U(a))/σ^2)

        # Return the escape rate
        return 1.0::Float64/(C*exp(-E))
end
