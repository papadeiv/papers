%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Fourier pseudo-spectral differentiation matrices - x =[0,L]
% Inputs: 2*L - length of domain, N - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Dx,D2x,w] = Compute_1D_Laplacian_fourier_neumann(N,L)

% Fourier mesh
M = 2*(N-1); dt = 2*pi/M; t = dt*(1:M)'; M2 = M/2;
column = [0 .5*(-1).^(1:M-1).*cot((1:M-1)*dt/2)];
Dt = toeplitz(column,column([1 M:-1:2]));
  
D2t = toeplitz([-pi^2/(3*dt^2)-1/6 ...
                 .5*(-1).^(2:M)./sin(dt*(1:M-1)/2).^2]);
  
% rewrite matrix for 0..pi reflect
semiDt = zeros(M2+1);
semiDt(:,1) = Dt(M2:M,M2);
semiDt(:,2:M2) = Dt(M2:M,M2-1:-1:1)+Dt(M2:M,M2+1:M-1);
semiDt(:,M2+1)=Dt(M2:M,M);
  
semiD2t = zeros(M2+1);
semiD2t(:,1) = D2t(M2:M,M2);
semiD2t(:,2:M2) = D2t(M2:M,M2-1:-1:1)+D2t(M2:M,M2+1:M-1);
semiD2t(:,M2+1)=D2t(M2:M,M);

% Rescale differentiation matrices to [-Lx,Lx]
semiDt  = (pi/L)*semiDt;        
semiD2t = (pi/L)^2*semiD2t;

x = L*(t - pi)/pi;
x = x(M2:M);
%% 1D
Dx = semiDt;
D2x= semiD2t;

ww = 2*L*ones(1,M)/M;   % trapezoidal weights for intergration int = w*u
semiw = zeros(1,M2+1);
semiw(1) = ww(M2); semiw(2:M2) = ww(M2+1:M-1) + ww(M2-1:-1:1);semiw(M2+1) = ww(M);
w = semiw;
