%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 4th order 1D finite-difference matrices - Neumann bcs on x=[0,L]
% Inputs: L - length of domain, N - number of mesh points
% Outputs: x - mesh, Dx - 1st differentiation matrix, D2x - 2nd order diff
% matrix, D4x - 4th order diff matrix, w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Dx,D2x,D4x,w] = Compute_1D_5ptLaplacian_finite_difference(N,L)

hx = L/(N-1);       % step size
x = (0:N-1)'*hx;    % x mesh

%% Finite-difference matrices
ex = ones(N,1);

Dx = spdiags([ex -8*ex 0*ex 8*ex -ex],-2:2, N, N);
Dx(1,:)   = 0; Dx(2,2)   = 1;
Dx(N,:) = 0; Dx(N-1,N-1) = -1;
Dx = Dx/(12*hx);

D2x = spdiags([-ex 16*ex -30.*ex 16*ex -ex], -2:2, N, N);
D2x(1,2)= 32; D2x(1,3)= -2;
D2x(2,1)= 16; D2x(2,2)= -31; 
D2x(N,:) = D2x(1,N:-1:1);
D2x(N-1,:) = D2x(2,N:-1:1);

D2x = D2x/(12*hx^2);

D4x = spdiags([-ex 12*ex -39*ex 56*ex -39*ex +12*ex -ex],-3:3, N, N);
D4x(1,2) = -78; D4x(1,3) = 24; D4x(1,4) = -2;
D4x(2,2) =  68; D4x(2,3) = -40;
D4x(3,2) = -40;
D4x(N,:) = D4x(1,N:-1:1);
D4x(N-1,:) = D4x(2,N:-1:1);
D4x(N-2,:) = D4x(3,N:-1:1);

D4x = D4x/(6*hx^4);

% w = [1, 2*ones(1,N-2), 1]*hx;   % trapezoidal weights for intergration int = w*u
w = [1, 2*ones(1,N-2)+2*mod([1:N-2],2),1]*hx/3; % simpson's rule weights for integration int = w*u

end
