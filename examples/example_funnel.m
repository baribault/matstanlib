
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script demonstrates "Neal's funnel" from the Stan manual.  
% this is useful for demonstrating matstanlib's diagnostic plots.  
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib 
%
% (c) beth baribault 2020 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_funnel';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 1000;     %how many iterations during warmup?
nIterations = 1000;     %how many iterations after warmup?

%an absolute path to a location for (temporary) Stan output files
% workingDir = '/my/custom/dir/';
workingDir = [pwd filesep 'wdir_' modelName];

%delete all temp files after stan finishes and selected output is saved?
%(this helps prevent compilation conflicts, etc.)
cleanUp = true;

%% prepare for output
fprintf('\n\n**********\n\npreparing to run the %s model.\n', ...
    upper(modelName))

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
stanfile = [workingDir modelName '.stan'];
stanfileID = fopen(stanfile,'w');
fprintf(stanfileID,'%s\n',modelCode{:});
fclose(stanfileID);

%% compile the model
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanfile);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%% run the model
tic
fprintf('\nrunning the model ... \n\n**********\n\n')
fit =  sm.sampling('file',  stanfile, ...
            'model_name',   modelName, ...
            'sample_file',  modelName, ...
            'chains',       nChains, ...
            'warmup',       nWarmup, ...
            'iter',         nIterations, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
stanSummary = fit.print;
fprintf('\n**********\n\ndone!\n')
fprintf('\nsampling took %.2f seconds.\n',toc)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('MATLABStan',fit);
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
disp(rtable)
%print a report about all diagnostics
interpretdiagnostics(diagnostics,rtable)

%trace plots
tracedensity(samples,diagnostics)

%% plots
%autocorrelation plots
plotautocorr(samples,'y')
plotautocorr(samples,'x')

plotlp(diagnostics,rtable)
jointdensity(samples,'x','y',diagnostics)
parcoordivergent(samples,diagnostics)
plotdivergences(diagnostics)
