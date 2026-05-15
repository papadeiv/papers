%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 1D Fourier pseudo-spectral diff marices - x=[-L,L]
% Inputs: 2*L - length of domain, N - number of mesh points
% Outputs: r - mesh, Dx - differentiation matrix, D2x - 2nd order diff mat
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Dx,D2x,w] = Compute_1D_Laplacian_fourier(N,L)

% Fourier differentiation matrix first order for y between -pi and pi
h = 2*pi/N;  t = h*(1:N); x = L*(t-pi)/pi;
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
Dt = toeplitz(column,column([1 N:-1:2]));
Dtt = toeplitz([-pi^2/(3*h^2)-1/6 .5*(-1).^(2:N)./sin(h*(1:N-1)/2).^2]);

%% 1D
Dx = (pi/L)*Dt;
D2x= (pi/L)^2*Dtt;

w = 2*L*ones(1,N)/N;   % trapezoidal weights for intergration int = w*u
