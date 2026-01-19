function [F,J] = LLE(uu ,p , mesh_params, dim)

% Mesh and operators parameters
N     = mesh_params.N;
r     = mesh_params.r;
Dx    = mesh_params.Dx;
Dxx   = mesh_params.Dxx;

% Scalar solution field u(r,t)
u_re = uu(1:N);                                                                 % Real part of u(r,t)
u_im = uu(1+N:2*N);                                                             % Imaginary part of u(r,t)

% Continuation parameters
zeta  = p(1);
f  = p(2);
nu  = p(3);

% Assemble the diffusion matrix according to the value of the homotopy
if dim == 2
    % Inverse radial mesh
    r(1) = 1;
    R = spdiags(1./r,0,N,N);
    % Discretised radial laplacian
    Lap = Dxx + nu*R*Dx;
    % Update the boundary terms according to the limit case 2d_rr at r = 0 
    Lap(1,1) = 2*Dxx(1,1); 
    Lap(1,2) = 2*Dxx(1,2);
else
    Lap = Dxx;
end

% Set the epsilon parameter to 1 (as per in Knobloch 2018)
epsilon = 1.0;

% Define the N-dimensional vector field
F_re = -Lap*u_re + zeta*u_re - (u_re.^2 + u_im.^2).*u_re + epsilon*u_im;        % Real part of f(u)
F_im = -Lap*u_im + zeta*u_im - (u_re.^2 + u_im.^2).*u_im + epsilon*(f - u_re);  % Imaginary part of f(u)

% Concatenate the vectors to return a 2N-dimensional vector field
F = [F_re; F_im];

% Compute the Jacobian around the supplied solution
if nargout > 1
    J11 = -Lap + spdiags(zeta*ones(N,1) - (u_re.^2 + u_im.^2) - (2*u_re.^2),0,N,N);
    J12 = spdiags(+epsilon*ones(N,1) - 2*u_im.*(u_re),0,N,N);
    J21 = spdiags(-epsilon*ones(N,1) - 2*u_re.*(u_im),0,N,N);
    J22 = -Lap + spdiags(zeta*ones(N,1) - (u_re.^2 + u_im.^2) - (2*u_im.^2),0,N,N);
    J = [J11 J12
         J21 J22];
end