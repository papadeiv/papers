% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function [V,LAMBDA] = ComputeSpectrum(u,p)

  %% Compute linear operators
  [~,J] = my_F(u,p);

  %% Call direct eigenvalue solver
  [V,LAMBDA] = eig(full(J));

end
