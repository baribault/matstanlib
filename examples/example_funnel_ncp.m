
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

% close all
clear
clc

%% inputs

modelName = 'example_funnel';

parameterization = 'centered';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 250;      %how many iterations during warmup?
nIterations = 2000;     %how many iterations after warmup?

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

%% data
disp('simulating data ... ')

nSubjects = 25;
nTrials = 10;
nData = nSubjects*nTrials;

Subject = repelem(1:nSubjects,nTrials)';

%true parameters
mu = 0;
sigma = 0.1;
nu = normrnd(mu,sigma,[nSubjects 1]);
omega = 1;

trueValues = collecttruevalues(mu,sigma,nu,omega);

%data
x = NaN(size(Subject));
for s = 1:nSubjects
    subjData = normrnd(nu(s),omega,[nTrials 1]);
    x(Subject==s) = subjData;
    
end

%collect data for Stan
dataStruct = struct('N',nData,'S',nSubjects,'Subject',Subject,'x',x);

disp('done!')

%% centered parameterization

%%% model specification %%%
modelCode = {
    'data { '
    '  int N;'
    '  int S;'
    '  int Subject[N];'
    '  real x[N];'
    '}'
    'parameters { '
    '  real mu;'
    '  real<lower=0> sigma;'
    '  real nu[S];'
    '  real<lower=0> omega;'
    '}'
    'model { '
    '  mu ~ normal(0,sqrt(10));'
    '  sigma ~ gamma(2,1);'
    '  '
    '  for (s in 1:S) '
    '    nu[s] ~ normal(mu,sigma);'
    '  omega ~ gamma(2,1);'
    '  '
    '  '
    '  for (n in 1:N) '
    '    x[n] ~ normal(nu[Subject[n]],omega);'
    '}'
};

%write the model code to a .stan file
stanfile = [workingDir modelName '.stan'];
stanfileID = fopen(stanfile,'w');
fprintf(stanfileID,'%s\n',modelCode{:});
fclose(stanfileID);

%%% compile the model %%%
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanfile);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%%% run the model %%%
tic
fprintf('\nrunning the model ... \n\n**********\n\n')
fit =  sm.sampling('file',  stanfile, ...
            'model_name',   modelName, ...
            'sample_file',  modelName, ...
            'data',         dataStruct, ...
            'chains',       nChains, ...
            'warmup',       nWarmup, ...
            'iter',         nIterations, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
stanSummary_c = fit.print;
fprintf('\n**********\n\ndone!\n')
fprintf('\nsampling took %.2f seconds.\n',toc)

%%% extract from the StanFit object %%%
[samples_c,diagnostics_c] = extractsamples('MATLABStan',fit);
%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); end

%%% diagnostics %%%
%calculate convergence diagnostics based on the posterior samples
rtable_c = rhattable(samples_c);
disp(rtable_c)
%print a report about all diagnostics
interpretdiagnostics(diagnostics_c,rtable_c)

%% non-centered parameterization

%%% model specification %%%
modelCode = {
    'data { '
    '  int N;'
    '  int S;'
    '  int Subject[N];'
    '  real x[N];'
    '}'
    'parameters { '
    '  real mu;'
    '  real<lower=0> sigma;'
    '  real nu_raw[S];'
    '  real<lower=0> omega;'
    '}'
    'transformed parameters { '
    '  real nu[S];'
    '  for (s in 1:S) '
    '    nu[s] = nu_raw[s]*omega + mu;'
    '}'
    'model { '
    '  mu ~ normal(0,sqrt(10));'
    '  sigma ~ gamma(2,1);'
    '  '
    '  for (s in 1:S) '
    '    nu_raw[s] ~ normal(0,1);'
    '  omega ~ gamma(2,1);'
    '  '
    '  '
    '  for (n in 1:N) '
    '    x[n] ~ normal(nu[Subject[n]],omega);'
    '}'
};

%write the model code to a .stan file
stanfile = [workingDir modelName '.stan'];
stanfileID = fopen(stanfile,'w');
fprintf(stanfileID,'%s\n',modelCode{:});
fclose(stanfileID);

%%% compile the model %%%
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanfile);
sm.compile();
fprintf('done!\n')
fprintf('compiling took %.2f seconds.\n',toc)

%%% run the model %%%
tic
fprintf('\nrunning the model ... \n\n**********\n\n')
fit =  sm.sampling('file',  stanfile, ...
            'model_name',   modelName, ...
            'sample_file',  modelName, ...
            'data',         dataStruct, ...
            'chains',       nChains, ...
            'warmup',       nWarmup, ...
            'iter',         nIterations, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
stanSummary_nc = fit.print;
fprintf('\n**********\n\ndone!\n')
fprintf('\nsampling took %.2f seconds.\n',toc)

%%% extract from the StanFit object %%%
[samples_nc,diagnostics_nc] = extractsamples('MATLABStan',fit);
%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%%% diagnostics %%%
%calculate convergence diagnostics based on the posterior samples
rtable_nc = rhattable(samples_nc);
disp(rtable_nc)
%print a report about all diagnostics
interpretdiagnostics(diagnostics_nc,rtable_nc)

%% plots

%centered parameterization
plotdivergences(diagnostics_c)
jointdensity(samples_c,'mu','sigma',diagnostics_c)
plotrecovery(getsamplestats(samples_c),trueValues)

%non-centered parameterization
plotdivergences(diagnostics_nc)
jointdensity(samples_nc,'mu','sigma',diagnostics_nc)
plotrecovery(getsamplestats(samples_nc),trueValues)
