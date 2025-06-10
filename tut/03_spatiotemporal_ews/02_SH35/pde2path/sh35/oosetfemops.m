function p=oosetfemops(p) % for SH as 2nd order system, hence singular p.mat.M  
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % scalar laplacian and mass 
p.mat.Dx=makeDx(p); % first order differentiation needed for H 
% Stiffness matrix
p.mat.K=[[0*K -K];[K M]];
% Mass matrix (singular)
p.mat.M=[[M 0*M];[0*M 0*M]];