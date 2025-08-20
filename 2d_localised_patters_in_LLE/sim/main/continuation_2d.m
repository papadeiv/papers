clear all, close all, clc;
addpath('./Secant_Continuation/');

L = 80;                                         % Domain
N = 500;                                        % Number of nodes
hx = L/(N-1);                                   % Discretisation size
dim = 2;                                        % Spatial dimension of the domain

% Radial mesh
r = (0:N-1)'*hx;

% Parameter vector
p(1) = 1.5;                                     % zeta == theta (fixed)
p(2) = 0.0;                                     % f == rho      (pattern continuation)
p(3) = 1.0;                                     % nu            (homotopy continuation)

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
sol = load("./cont_nu_step_020/solution_0001320.mat");                % Import the struct from the mat file
u0 = sol.u;                                                           % Import solution from the 1D continuation
u0_re = u0(1:N);                                                      % Extract the real part
u0_im = u0(N+1:2*N);                                                  % Extract the imaginary part
plot(r,u0_re,'b',r,u0_im,'r');  drawnow;

% Extract the current value of rho out of the initial guess
p(2) = sol.p(2);

% Solve the 0-problem to find the initial guess
options = optimset('Display','iter','DerivativeCheck','off','Jacobian','off','MaxIter',100,'Algorithm','levenberg-marquardt');
[u_converged,fval,exitflag,output,jacobian] = fsolve(@(u) LLE(u,p,mesh_params,dim), u0, options);

% Plot the converged initial guess
figure;plot(r,u_converged(1:N));drawnow;

%% Numerical continuation in rho using Daniele's code
% Various function handles
problemHandle            = @(u,p) LLE(u,p,mesh_params,dim);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p,mesh_params);
computeEigenvaluesHandle = @(u,p) ComputeSpectrum(u,p,mesh_params,dim);
plotSpetcrumHandle       = [];
stepperPars.iContPar      = 2;     % continuation parameter index p(2)
stepperPars.s0            = -1e-3; % initial arclength step size - minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;  % smallest arclength step size allowed
stepperPars.sMax          = 1e-2;  % largest arclength step size allowed
stepperPars.pMin          = 0.0;   % minimum paramter value
stepperPars.pMax          = 4.0;   % maximum parameter value
stepperPars.maxSteps      = 20000; % maximum number of continuation steps
stepperPars.nPrint        = 1;     % print solution stats every nPrint steps
stepperPars.nSaveSol      = 10;    % save solution every nSaveSol steps
stepperPars.finDiffEps    = 1e-10; % finite diff step for arclength continuation
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 2e1;   % number of nonlinear iterations 
stepperPars.dataFolder    = 'cont_rho_from_nu_step_020_3'; % folder to save data
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1; % plot Id from the SolutionMeasures.m
stepperPars.uzstop= @(v1,v0,val) 0; % stopping criterion
branch = SecantContinuation(problemHandle,u_converged,p,stepperPars);

%% Plot the radially symmetric solution in 2d

% Get the total number of continuation steps
sol_branch = load("./cont_rho_from_nu_step_020_3/branch.mat");
N_steps = size(sol_branch.branch,1)-1;

% Get the solution at the end of the continuation
filename = sprintf('./cont_rho_from_nu_step_020_3/solution_000%04d.mat', N_steps);
solution = load(filename);
u = solution.u;

% Plot the solution in 2D
PlotSolution_2d(u, mesh_params)