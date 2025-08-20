function F = SolutionMeasures(step, uu, p, mesh_params)

N = mesh_params.N;
u_re = uu(1:N);
u_im = uu(N+1:2*N);

% Compute the homogeneous steady-state
[U_0, V_0] = ComputeHomogeneousSteadyState(p(1), p(2));

% Compute the norm of the localised solution without the background state
L2norm = sqrt(mean((u_re - U_0).^2 + (u_im - V_0).^2));
F = [L2norm];

end
