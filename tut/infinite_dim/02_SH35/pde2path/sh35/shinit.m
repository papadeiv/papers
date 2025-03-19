function p=shinit(p,nx,lx,ly,ndim,par,varargin) % Swift-Hohenberg as 2 component system, with 
% singular M=diag(M,0) mass-matrix to compute correct eigenvalues 

% Check the fields of the struct
% file:///home/dpap666/Libraries/MATLAB/pde2path/html/index/stanparam.html

% Standard problem init
p=stanparam(p);
% Number of equations
p.nc.neq=2;
% Number of spatial dimensions
p.ndim=ndim; 
% Number of eigenvalues of the Jacobian of G(u)
p.nc.neig=20; 
% Function handle (fuha) of the non-linear LHS G(u)
p.fuha.sG=@sG; 
% Function handle (fuha) of the Jacobian of G(u)
p.fuha.sGjac=@sGjac; 
% Don't calculate stability EVals
p.sw.spcalc=0;
% Component of branch for plotting
p.plot.bpcmp=0; 
% Stepsize of the continuation starting from p.u in direction p.tau
p.sol.ds=0.001; 
% Max stepsize
p.sol.dsmax=0.01; 
% Residual accuracy
p.pm.resfac=1e-1; 
% Check for bifurcations using spcalc field
p.sw.bifcheck=2; 
% Normal continuation
p.sw.spcont=0;
% Do calculate stability EVals
p.sw.spcalc=1;
% Switch for simple Jacobian and residuals
p.sw.sfem=-1; 
% Fold detection on
p.sw.foldcheck=1; 

% Generate PDE object
pde=stanpdeo1D(lx,2*lx/nx); 
% Spatial domain (Om=Ω)
p.Om=2*lx;
% Branch output
p.fuha.outfu=@shbra1d;

%
p.np=pde.grid.nPoints;
%
p.pdeo=pde;
%
p.nu=p.np*p.nc.neq;
%
p.sol.xi=1/p.nu;
% Parameter range for the bifurcation diagram
p.nc.lammin=-1.25; 
p.nc.lammax=1.25;
%
p.nc.ds=0.001; 
%
p.nc.dsmax=0.01; 
u=0*ones(p.np,1); 
v=u; 
u0=[u v]; 
p.u=u0(:); 
p.u=[p.u; par]; 
p.nc.ilam=1; 
p.file.smod=5; 
%screenlayout(p); 
p=setfemops(p); 
p.nc.nsteps=3000;