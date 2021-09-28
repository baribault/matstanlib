
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%     EXAMPLE_FUNNEL    %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script showcases the "funnel" geometry characteristic of
% hierarchical models. 
% 
% the model specifications are borrowed from Betancourt & Girolami (2016).
% 
% first, a simple linear model is implemented using a centered 
% parameterization.  in this version, the bottom of the funnel is unable 
% to be explored fully.  the problematic sampling is (typically) flagged by 
% the automated diagnostics.  
% next, the same model is implemented using a non-centered parameterization.
% in this version, the funnel is able to be explored fully, and the sampling
% is of good quality.  
% 
% this example emphasizes how matstanlib's automated diagnostics can help to 
% detect problems with sampling, and how the diagnostic plotting functions 
% help to identify the root of the issue.  
% it also emphasizes how one can confirm that a solution (here, the use of 
% reparameterization) remedies the problem by observing its effects on a 
% second set of the same diagnostics and plots.
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
%
% Reference:    Betancourt, M. & Girolami, M. (2016).  Hamiltonian Monte 
%                   Carlo for hierarchical models.  In S. Upadhyay, A.
%                   Loganathan & U. Singh (Eds.), Current Trends in 
%                   Bayesian Methodology with Applications (pp. 79–-101).
%                   Chapman & Hall/CRC.
% 
% (c) beth baribault 2020 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_funnel';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 500;      %how many iterations during warmup?
nIterations = 1500;     %how many iterations after warmup?

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
fprintf('\nsimulating data ... ')

nSubjects = 20;
nTrials = 15;
nData = nSubjects*nTrials;

Subject = repelem(1:nSubjects,nTrials)';

%true parameters
mu = 0;
sigma = 0.1;
nu = normrnd(mu,sigma,[nSubjects 1]);
omega = 1;

trueValues = collecttruevalues(mu,sigma,nu);

%data
x = NaN(size(Subject));
x_sd = NaN([1 nSubjects]);
for s = 1:nSubjects
    subjData = normrnd(nu(s),omega,[nTrials 1]);
    x(Subject==s) = subjData;
    x_sd(s) = std(subjData);
end

%collect data for Stan
dataStruct = struct('N',nData,'S',nSubjects,'Subject',Subject,'x',x,'x_sd',x_sd);

disp('done!')

%% model specifications

%%% centered parameterization %%%
modelCode_centered = {
    'data { '
    '  int N;'
    '  int S;'
    '  int Subject[N];'
    '  real x[N];'
    '  real x_sd[S];'
    '}'
    'parameters { '
    '  real mu;'
    '  real<lower=0> sigma;'
    '  real nu[S];'
    '}'
    'model { '
    '  mu ~ normal(0,sqrt(10));'
    '  sigma ~ gamma(2,1);'
    '  '
    '  for (s in 1:S) '
    '    nu[s] ~ normal(mu,sigma);'
    '  '
    '  for (n in 1:N) '
    '    x[n] ~ normal(nu[Subject[n]],x_sd[Subject[n]]);'
    '}'
};

%%% non-centered parameterization %%%
modelCode_noncentered = {
    'data { '
    '  int N;'
    '  int S;'
    '  int Subject[N];'
    '  real x[N];'
    '  real x_sd[S];'
    '}'
    'parameters { '
    '  real mu;'
    '  real<lower=0> sigma;'
    '  real nu_raw[S];'
    '}'
    'transformed parameters { '
    '  real nu[S];'
    '  for (s in 1:S) '
    '    nu[s] = nu_raw[s]*sigma + mu;'
    '}'
    'model { '
    '  mu ~ normal(0,sqrt(10));'
    '  sigma ~ gamma(2,1);'
    '  '
    '  for (s in 1:S) '
    '    nu_raw[s] ~ normal(0,1);'
    '  '
    '  for (n in 1:N) '
    '    x[n] ~ normal(nu[Subject[n]],x_sd[Subject[n]]);'
    '}'
};

%% fit both versions of the model
for m = 1:2
    if m==1
        modelCode = modelCode_centered;
    elseif m==2
        modelCode = modelCode_noncentered;
    end
    %write the model code to a .stan file
    stanFilePath = writestanfile(modelCode,modelName,workingDir);
    
    %compile the model
    tic
    fprintf('\ncompiling the model ... ')
    sm = StanModel('file',stanFilePath);
    sm.compile();
    fprintf('done!\n')
    fprintf('compiling took %.2f seconds.\n',toc)
    
    %run the model
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
    
    %extract samples, run diagnostics
    if m==1
        %extract from the StanFit object
        [samples_c,diagnostics_c] = extractsamples('MATLABStan',fit);
        %clean up after MATLABStan
        clearvars fit
        %clean up after Stan
        if cleanUp, delete([workingDir '*']); end

        %%% diagnostics %%%
        %compute posterior sample-based diagnostics and summary statistics
        posteriorTable_c = mcmctable(samples_c);
        disp(posteriorTable_c)
        %print a report about all MCMC diagnostics
        interpretdiagnostics(diagnostics_c,posteriorTable_c)
    elseif m==2
        %extract from the StanFit object
        [samples_nc,diagnostics_nc] = extractsamples('MATLABStan',fit);
        %clean up after MATLABStan
        clearvars fit
        %clean up after Stan
        if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

        %%% diagnostics %%%
        %compute posterior sample-based diagnostics and summary statistics
        posteriorTable_nc = mcmctable(samples_nc);
        disp(posteriorTable_nc)
        %print a report about all MCMC diagnostics
        interpretdiagnostics(diagnostics_nc,posteriorTable_nc)
    end
end

%% plots

%the centered parameterization generally leads to divergences ...
plotdivergences(diagnostics_c)
%... which concentrate at the bottom of the "funnel" 
%    (where mu & sigma are both small)
jointdensity(samples_c,'mu','sigma',diagnostics_c)

%we can resolve this by using a non-centered parameterization instead, 
%which does not lead to divergences ...
plotdivergences(diagnostics_nc)
% ... because the sampler is no longer prevented reaching the bottom of 
%     the funnel, and so can fully explore the posterior
jointdensity(samples_nc,'mu','sigma',diagnostics_nc)

