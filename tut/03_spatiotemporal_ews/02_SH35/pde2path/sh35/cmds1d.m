%% Clear the results of the previous run
keep pphome; close all; clc;
problem=[]; 
%% Initialise the model and specify the numerical settings
% Endpoints of the spatial domain
lx=32*pi;
% Number of nodes of the discretised grid
nx=round(100*lx);
% Fictitious values for the y-variable (which does not exist in 1D)
ly=0.1;
ny=1;
% Number of spatial dimensions
ndim=1;
% Initial value of the bifurcation parameter
lam=-1.0;
% Initial value of the non-bifurcating parameter (kept fixed)
nu=2;
% Parameter vector
par=[lam; nu];
% Initialise the problem struct
problem=shinit(problem,nx,lx,ly,ndim,par);
huclean(problem);
% Create the output directory
problem=setfn(problem,'1D/run2');
% Look for n bifurcations in the parameter range
n = 100;
problem=findbif(problem,n);
%% First continuation run at lambda = 0 (periodic solution)
problem=swibra('1D/run0','bpt1','1D/run1',0.001);
problem.sol.ds=0.01; 
problem.pm.resfac=1e-1; 
problem=cont(problem,1000);
%% First continuation run at lambda = 1 (homogeneous solutions)
problem=swibra('1D/run2','bpt55','1D/run3',0.001);
problem.sol.ds=0.01; 
problem.pm.resfac=1e-1; 
problem=cont(problem,1000);
%% Plot the bifurcation diagram
pcmp=3; 
figure(3); 
clf; 
% Ground-state branch
plotbra('1D/run0','pt625','ms',0,'tyst','-','lwst',3,'lwun',1);
% Homogeneous solutions' branch
plotbra('1D/run3','pt270','ms',0,'tyst','-','lwst',3,'lwun',1);
% Spatially-periodic solutions' branch
plotbra('1D/run1','pt393','ms',1,'tyst','-','lwst',3,'lwun',1,'cl','b'); 
xlabel('\mu'); 
ylabel('||u||_2');
