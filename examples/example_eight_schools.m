
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE_EIGHT_SCHOOLS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script implements the classic "eight schools" example from the 
% Bayesian Data Analysis textbooks (with a non-centered parameterization). 
% 
% this familiar model specification is a useful context you to test out 
% matstanlib's capabilites.  
% 
% this example demonstrates how matstanlib is not only useful for diagnostic 
% purposes, but also for supporting a good Bayesian workflow that emphasizes 
% visual exploration of your model output.  
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
% 
% Reference:    Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., 
%                   Vehtari, A., & Rubin, D.B. (2013). Bayesian data 
%                   analysis (3rd ed.). Chapman & Hall/CRC.  
% 
% (c) beth baribault 2021 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_eight_schools';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 1000;     %how many iterations during warmup?
nIterations = 1000;     %how many iterations after warmup?

%an absolute path to a (temporary) location for Stan output files:
workingDir = [pwd filesep 'wdir_' modelName];

%delete temporary files after Stan finishes and output has been returned?
%(this helps prevent compilation conflicts, etc.)
cleanUp = true;

%% setup
fprintf('\n\n************\n\npreparing to run the %s model.\n',modelName)

%working directory
if ~isequal(workingDir(end),filesep), workingDir(end+1) = filesep; end
fprintf(['\nStan files will be saved to the working directory:\n' ...
    '> %s\n'], workingDir)
if ~exist(workingDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(workingDir)
    fprintf('done.\n')
end

%% data
fprintf('\nloading data ... ')

%collect data for Stan
dataStruct = struct('J',8, ...
                    'y',[28 8 -3 7 -1 1 18 12],...
                    'sigma',[15 10 16 11 9 11 10 18]);

fprintf('done!\n')

%% model specification

modelCode = {
    '//from Section 5.5 of Gelman et al (2003)'
    'data { '
    '  int<lower=0> J;              //number of schools'
    '  real y[J];                   //estimated treatment effect'
    '  real<lower=0> sigma[J];      //standard error of the effect estimate'
    '}'
    'parameters { '
    '  real mu;                     //population treatment effect'
    '  real<lower=0> tau;           //standard deviation in treatment effects'
    '  vector[J] eta;               //unscaled deviation from mu by school'
    '}'
    'transformed parameters {'
    '  vector[J] theta;             //school treatment effects'
    '  for (j in 1:J) '
    '    theta[j] = mu + tau * eta[j];'
    '}'
    'model { '
    '  eta ~ normal(0,1);'
    '  y ~ normal(theta,sigma);'
    '}'
};

%write the model code to a .stan file
stanFilePath = writestanfile(modelCode,modelName,workingDir);

%% compile the model
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanFilePath);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%% run the model
tic
fprintf('\nrunning the model ... \n\n************\n\n')
fit =  sm.sampling('file',          stanFilePath, ...
                   'model_name',    modelName, ...
                   'sample_file',   modelName, ...
                   'verbose',       true, ...
                   'chains',        nChains, ...
                   'warmup',        nWarmup, ...
                   'iter',          nIterations, ...
                   'data',          dataStruct, ...
                   'working_dir',   workingDir);
fit.block();
[stanSummaryTxt,stanSummary] = fit.print('sig_figs',5);
fprintf('\n************\n\ndone!\n')
fprintf('\nsampling took %.2f seconds.\n',toc)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('MATLABStan',fit);
parameters = fieldnames(samples);
instances = getparaminstances([],samples);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%% diagnostic reports & plots
%compute posterior sample-based diagnostics and summary statistics
posteriorTable = mcmctable(samples);
%print a report about all MCMC diagnostics
interpretdiagnostics(diagnostics,posteriorTable)

%trace plots/rank plots
rankplots(samples,{'mu','tau'})

%% 
