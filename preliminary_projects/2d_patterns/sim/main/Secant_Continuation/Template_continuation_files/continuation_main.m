% continuation template files

clear all, close all, clc; % clear workspace, clear figures, clear screen

addpath('../Secant_Continuation/'); % this will need changing depending on where you put this in a folder

% define p a vector
% define initial guess u0

%% Continue in the bifurcation parameter Daniele's continuation code
% Various function handles
problemHandle            = @(u,p) my_F(u,p);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle);
branchVariablesHandle    = @(step,u,p) SolutionMeasures(step,u,p);
computeEigenvaluesHandle = @(u,p) ComputeSpectrum(u,p);
plotSpetcrumHandle       = @(lambda,p,parent) PlotSpectrum(lambda,p,parent);
stepperPars.iContPar      = 1;     % continuation parameter index p(1)
stepperPars.s0            = 0.01; % initial arclength step size - minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;  % smallest arclength step size allowed
stepperPars.sMax          = .1;    % largest arclength step size allowed
stepperPars.pMin          = -1.0;  % minimum paramter value
stepperPars.pMax          = 2;     % maximum parameter value
stepperPars.maxSteps      = 20000; % maximum number of continuation steps
stepperPars.nPrint        = 1;     % print solution stats every nPrint steps
stepperPars.nSaveSol      = 100;   % save solution every nSaveSol steps
stepperPars.finDiffEps    = 1e-7;  % finite diff step for arclength continuation
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);
stepperPars.optNonlinIter = 10;     % number of nonlinear iterations 
stepperPars.dataFolder    = 'Run_1'; % folder to save data
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 1; % plot Id from the SolutionMeasures.m
stepperPars.uzstop= @(v1,v0,val) uzstop_val(v1,v0,1); % stopping criterion
branch = SecantContinuation(problemHandle,u0,p,stepperPars);
