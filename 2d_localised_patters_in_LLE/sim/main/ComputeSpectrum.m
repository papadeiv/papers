% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function [V,LAMBDA] = ComputeSpectrum(uu,p,mesh_params,dim)

  % Compute the linear operators and the Jacobian
  [~,J] = LLE(uu,p,mesh_params,dim);
  [M, N] = size(J);

  % Because of complex system, we need to modify the Jacobian
  temp = J;
  J(1:M/2,:) = J(M/2+1:end,:);
  J(M/2+1:end,:) = temp(1:M/2,:);
  J(M/2+1:end,:) = -J(M/2+1:end,:);

  % Get the spectrum of the modified Jacobian
  [V,LAMBDA] = eig(full(J));
end
