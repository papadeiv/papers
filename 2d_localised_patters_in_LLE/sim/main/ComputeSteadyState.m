function [U,V] = ComputeSteadyState(theta, epsilon, wavenumber, phase, x)

% U = real part of u(r,t)
% V = imaginary part of u(r,t)

% Compute rho_c where the HH bifurcation occurs
rho_c = sqrt(1 + (1 - theta)^2);
% Compute the perturbation from rho_c
delta_c = ((theta - 2)^2)/(2*rho_c);
% Compute the neighbourhood of rho_c
rho = rho_c - (epsilon^2)*delta_c;

% Compute the HSS
U_c = rho_c/(1 + (1-theta)^2);
V_c = ((1-theta)*rho_c)/(1 + (1-theta)^2);

% Compute the leading-order perturbation
prefactor_1 = (rho - rho_c)/((theta^2 - 2*theta + 2)*(theta - 2));
U_2 = prefactor_1*(theta^2);
V_2 = prefactor_1*(-theta^2 -theta + 2);

% Compute the state-dependent correction
C_1 = -(2*(theta^2 - 2*theta + 2))/(theta - 2);
C_2 = (2*((theta^2 - 2*theta + 2)^(3/2)))/((theta - 2)^4);
C_3 = (4*((theta^2 - 2*theta + 2)^2)*(30*theta - 41))/(9*((theta - 2)^6));
a = theta/(2-theta);
f = sech(sqrt(C_2*(rho_c - rho)/(C_1))*x).*cos(wavenumber*x + phase);
U1 = 2*a*sqrt(2*(C_2*(rho_c - rho))/(C_3))*f;
V1 = 2*sqrt(2*(C_2*(rho_c - rho))/(C_3))*f;

% Assemble the initial guess
U = U_c + U1 + U_2;
V = V_c + V1 + V_2;

end