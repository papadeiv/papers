function [U,V] = ComputeHomogeneousSteadyState(theta, rho)

% U = real part of u(r,t)
% V = imaginary part of u(r,t)

% Compute the squared amplitude
I_0 = 1.0;

% Compute the HSS
U = rho/(1 + (I_0-theta)^2);
V = ((I_0-theta)*rho)/(1 + (I_0-theta)^2);

end