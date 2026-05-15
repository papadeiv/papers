clear all, close all, clc;
addpath('./Secant_Continuation/');

L = 80;                                         % Domain
N = 500;                                        % Number of nodes
hx = L/(N-1);                                   % Discretisation size
dim = 1;                                        % Spatial dimension of the domain

% Radial mesh
r = (0:N-1)'*hx;
% Inverse radial mesh

% PDE parameters
theta = 1.5;
rho_c = sqrt(1 + (1 - theta)^2);

% Perturbation parameters
delta_c = (theta - 2)^2/2*rho_c;
rho_0 =  1.11753;
epsilon = sqrt((rho_c-rho_0)/delta_c);

% Parameter vector
p(1) = 1.5;                                     % zeta == theta (fixed)
p(2) = rho_c - (epsilon^2)*delta_c;             % f == rho      (pattern continuation)
p(3) = 0.0;                                     % nu            (homotopy continuation)

%% Discretise the PDE (using finite differences) and set-up the 0-problem
ex = ones(N,1);

% Discretised gradient
Dx = spdiags([-ex 0*ex ex],-1:1, N, N);
Dx(1,2) = 0; Dx(N,N-1) = 0;                     % Impose Neumann bcs
Dx = Dx/(2*hx);

% Discretised cartesian Laplacian
Dxx = spdiags([ex -2*ex ex], -1:1, N, N);
Dxx(1,2)=2;Dxx(N,N-1)=2;                        % Impose Neumann bcs
Dxx = Dxx/(hx^2);

% Mesh and operator parameters
mesh_params.N = N; 
mesh_params.r = r;
mesh_params.Dx= Dx;
mesh_params.Dxx = Dxx;

% Get a good starting solution for the initial guess
phi = 0.0;                                                            % Phase (either 0 or pi)
k = sqrt(0.5);                                                        % Wavenumber (-Im(lambda))
[u0_re, u0_im] = ComputeSteadyState(p(1), epsilon, k, phi, r);        % Compute the non-homogenous steady-state (as provided by (27))
u0 = [u0_re(:); u0_im(:)];                                            % Combine real and imaginary parts into a single vector
plot(r,u0_re,'b',r,u0_im,'r');  drawnow;

% Solve the 0-problem to find the initial guess
options = optimset('Display','iter','DerivativeCheck','on','Jacobian','on','MaxIter',100,'Algorithm','levenberg-marquardt');
[u_converged,fval,exitflag,output,jacobian] = fsolve(@(u) LLE(u,p,mesh_params,dim), u0, options);

% Plot the converged initial guess
figure;plot(r,u_converged(1:N));drawnow;

%% Numerical continuation in rho using Daniele's code
% Various function handles
problemHandle            = @(u,p) LLE(u,p,mesh_params,dim);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p,mesh_params);
computeEigenvaluesHandle = @(u,p) ComputeSpectrum(u,p,mesh_params,dim);
plotSpetcrumHandle       = @(lambda,p,parent) PlotSpectrum(lambda,p,parent);
stepperPars.iContPar      = 2;     % continuation parameter index p(2)
stepperPars.s0            = -1e-3; % initial arclength step size - minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;  % smallest arclength step size allowed
stepperPars.sMax          = 5e-2;  % largest arclength step size allowed
stepperPars.pMin          = 1.0;   % minimum paramter value
stepperPars.pMax          = 1.5;   % maximum parameter value
stepperPars.maxSteps      = 20000; % maximum number of continuation steps
stepperPars.nPrint        = 1;     % print solution stats every nPrint steps
stepperPars.nSaveSol      = 10;    % save solution every nSaveSol steps
stepperPars.finDiffEps    = 1e-10; % finite diff step for arclength continuation
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','on',...
                                     'Jacobian','on',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 2e1;   % number of nonlinear iterations 
stepperPars.dataFolder    = 'cont_rho_0'; % folder to save data
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1; % plot Id from the SolutionMeasures.m
stepperPars.uzstop= @(v1,v0,val) 0; % stopping criterion
branch = SecantContinuation(problemHandle,u_converged,p,stepperPars);