
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this file is a skeleton of a script for running your own Stan      %%%
%%% models with Stan, MATLABStan, and matstanlib.                      %%%
%%% (c) beth baribault 2019 ---                           > matstanlib %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% a bare bones outline of a model fitting script
% 
% ...
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
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
nIterations = 1000;     %how many iterations after warmup?

%an absolute path to a location for (temporary) Stan output files
% workingDir = '\my\custom\folder\';    %%% windows
% workingDir = '/my/custom/folder/';    %%% macOS, linux
workingDir = [pwd filesep 'wdir_' modelName];

%delete all temp files after stan finishes and selected output is saved?
%(this helps prevent compilation conflicts, etc.)
cleanUp = true;

%% prepare for output
fprintf('\n\n**********\n\npreparing to run the ''%s'' model.\n', ...
    modelName)

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
% disp('simulating data ... ')
% disp('loading data ... ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% either simulate or load data here %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%collect data for Stan
dataStruct = struct('%%%%%%%%%%%');

% disp('done!')

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
    '  //predicitive distributions and other tracked quantites go here.'
    '  //...'
    '}'
};

%write the model code to a .stan file
stanFile = [workingDir modelName '.stan'];
stanFileID = fopen(stanFile,'w');
fprintf(stanFileID,'%s\n',modelCode{:});
fclose(stanFileID);
%copy the .stan file to current directory
copyfile(stanFile,pwd)

%% compile the model
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanFile);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%% run the model
tic
fprintf('\nrunning the model ... \n\n**********\n\n')
fit =  sm.sampling('file',  stanFile, ...
            'model_name',   modelName, ...
            'sample_file',  modelName, ...
            'data',         dataStruct, ...
            'chains',       nChains, ...
            'warmup',       nWarmup, ...
            'iter',         nIterations, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
stan_summary = fit.print;
fprintf('\n**********\n\ndone!\n')

runtime = datevec(seconds(toc));
fprintf(['sampling took %i days, %i hours, %i minutes, ' ...
    'and %.2f seconds.\n'],runtime(3:end)) %optimistically assuming your 
                                           %model takes < 1 month to run :)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('matlabstan',fit);
parameters = fieldnames(samples);
instances = getparaminstances([],samples);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp
    delete([workingDir '*'])
    rmdir(workingDir)
end

%% diagnostics
%calculate convergence diagnostics based on the posterior samples
rtable = rhattable(samples);
%print a report about all diagnostics
interpretdiagnostics(diagnostics,rtable)

%trace plots
tracedensity(samples,diagnostics)

%% parameter estimates & model comparison & other statistics
estimatedValues = getsamplestats(samples);
% estimatedValues = getsamplestats(samples,trueValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   compute other statistics here   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save output
save([pwd filesep modelName '.mat'])

%% generate plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   call plotting functions, etc.   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
