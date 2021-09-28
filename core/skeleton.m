
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this file is a skeleton of a script that you may use as a starting  %%%
%%% point for your own code with Stan, MATLABStan, and matstanlib.      %%%
%%%                                                                     %%%
%%% (c) beth baribault 2019 ---                            > matstanlib %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% a bare bones outline of a Bayesian model fitting script
% 
% ...
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/users/interfaces/cmdstan.html
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
%
% (this script was built from matstanlib's skeleton.m file.)

close all
clear
clc

%% inputs

modelName = 'my_model_name';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 1000;     %how many iterations during warmup?
nIterations = 2000;     %how many iterations after warmup?

%an absolute path to a location for saved model output, figures, etc.:
% outputDir = 'X:\\my\custom\dir\'; %%% windows
% outputDir = '~/my/custom/dir/';   %%% macOS, linux
outputDir = [pwd filesep 'stan_output_' modelName];

%an absolute path to a (temporary) location for Stan output files:
workingDir = [pwd filesep 'wdir_' modelName];

%delete temporary files after Stan finishes and output has been returned?
%(this helps prevent compilation conflicts, etc.)
cleanUp = true;

%% setup
fprintf('\n\n************\n\npreparing to run the %s model.\n',modelName)

%ensure all directories inputs end in a file separator
if ~isequal(outputDir(end),filesep), outputDir(end+1) = filesep; end
if ~isequal(workingDir(end),filesep), workingDir(end+1) = filesep; end

%output directory
fprintf(['\nmodel output and figures will be saved to the output directory:\n' ...
    '> %s\n'], outputDir)
if ~exist(outputDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(outputDir)
    fprintf('done.\n')
end

%working directory
fprintf(['\nStan files will be saved to the working directory:\n' ...
    '> %s\n'], workingDir)
if ~exist(workingDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(workingDir)
    fprintf('done.\n')
end

%% data
% disp('\nsimulating data ... ')
% disp('\nloading data ... ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% either simulate or load data here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%collect data for Stan
dataStruct = struct('%%%%%%%%%%%');

% disp('done!\n')

%% model specification

modelCode = {
    'data { '
    '  //input data variables are declared here.'
    '  //...'
    '}'
    ''
    'transformed data { '
    '  //transformations of data and fixed variables are declared here.'
    '  //...'
    '}'
    ''
    'parameters { '
    '  //parameters are declared here.'
    '  //...'
    '}'
    ''
    'transformed parameters {'
    '  //transformations and derived parameter declarations go here.'
    '  //...'
    '}'
    ''
    'model { '
    '  //priors and likelihood go here.'
    '  //...'
    '}'
    ''
    'generated quantities {'
    '  //predictive distributions and other quantites to track go here.'
    '  //...'
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
fprintf('\n************\n\ndone!\n\n')

runtime = datevec(seconds(toc));
fprintf(['sampling took %i days, %i hours, %i minutes, ' ...
    'and %.2f seconds.\n'],runtime(3:end)) %optimistically assuming your 
                                           %model takes < 1 month to run :)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('MATLABStan',fit);
parameters = fieldnames(samples);
instances = getparaminstances([],samples);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%% save output
save([pwd filesep modelName '.mat'])

%% diagnostic reports & plots
%compute posterior sample-based diagnostics and summary statistics
posteriorTable = mcmctable(samples);
%print a report about all MCMC diagnostics
interpretdiagnostics(diagnostics,posteriorTable)

%trace plots/rank plots
tracedensity(samples,'???',diagnostics)
rankplots(samples,'???')

%% parameter estimates & model comparison & other statistics
estimatedValues = getsamplestats(samples);
% estimatedValues = getsamplestats(samples,trueValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   call plotting functions, etc.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% other analyses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   compute other statistics here   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   call plotting functions, etc.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
