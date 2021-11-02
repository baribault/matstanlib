
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  EXAMPLE_FUNNEL_NEAL  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script implements "Neal's funnel" from the Stan manual. 
% 
% this example demonstrates some of matstanlib's plots that are useful 
% for diagnosing divergences.  
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/users/interfaces/cmdstan.html
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib 
%
% Reference: 
%     https://mc-stan.org/docs/2_27/stan-users-guide/reparameterization-section.html
% 
% (c) beth baribault 2020 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_neals_funnel';

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

%% model specification

modelCode = {
    'parameters { '
    '  real y;'
    '  real x;'
%     '  vector[3] x;'
    '}'
    'model { '
    '  y ~ normal(0,3);'
    '  x ~ normal(0,exp(y/2));'
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
tracedensity(samples,{},diagnostics)
rankplots(samples)

%% diagnostic plots
%autocorrelation plots
plotautocorr(samples,'y')
plotautocorr(samples,'x')

%bird's eye view
plotlp(diagnostics,posteriorTable)

%sure do have a lot of divergences there
plotdivergences(diagnostics)

%identifying the issue
jointdensity(samples,'x','y',diagnostics)
parcoordivergent(samples,diagnostics)
