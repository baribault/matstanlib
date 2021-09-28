
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  EXAMPLE_CORRELATION  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script implements a Pearson correlation as a Bayesian model. 
% 
% the model specification is an adaptation for Stan of the 'Correlation_1'
% model from Chapter 5.1 in Lee & Wagenmakers (2014; pp. 60--61). 
% 
% this example demonstrates good structure for MATLAB code that implements
% Bayesian models of any flavor: 
%           (1) a section just for inputs
%           (2) basic setup (e.g., setting paths for Stan output)
%           (3) loading/simulating data & data formatting
%           (4) model specification
%           (5) compiling and running the model via Stan
%           (6) extracting samples from the interface & cleaning up
%           (7) running diagnostic checks
%           (-) ... and so on. 
% it also emphasizes how matstanlib's plots support assessment of the 
% quality of parameter recovery in addition to a Bayesian hypothesis test.
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
% 
% Reference:    Lee, M.D. & Wagenmakers, E.J. (2014). Bayesian cognitive
%                   modeling: A practical course. Cambridge University 
%                   Press. 
% 
% (c) beth baribault 2019 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_correlation';

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
fprintf('\nsimulating data ... ')

%known features of the experiment
nObservations = 100;

%true parameter values
rho     = rand*2 - 1; %pearson correlation is [-1 1]
mu      = [3 5];
sigma   = [0.75 1.25];
E       = [sigma(1)^2             rho*sigma(2)*sigma(1); ...
           rho*sigma(1)*sigma(2)  sigma(2)^2          ];

%simulate observed data
x = mvnrnd(mu,E,nObservations);

%create a structure of true values
trueValues = collecttruevalues(mu,sigma,rho);

%collect data for Stan
data = struct('x',x,'N',nObservations);

fprintf('done!\n')

%% model specification

modelCode = {
    'data { '
    '  int N;                             //number of observations'
    '  vector[2] x[N];                    //"observed" data'
    '}'
    'parameters { '
    '  vector[2] mu;'
    '  vector<lower=0>[2] sigma;'
    '  real<lower=-1,upper=1> rho;        //pearson correlation'
    '}'
    'transformed parameters {'
    '  vector<lower=0>[2] tau;'
    '  cov_matrix[2] E;'
    '  tau[1] = pow(sigma[1],-2);         //sd vs. precision'
    '  tau[2] = pow(sigma[2],-2);         //sd vs. precision'
    '  E[1,1] = sigma[1]^2;               //variance   x1x1'
    '  E[1,2] = rho*sigma[2]*sigma[1];    //covariance x2x1'
    '  E[2,1] = rho*sigma[1]*sigma[2];    //covariance x1x2'
    '  E[2,2] = sigma[2]^2;               //variance   x2x2'
    '}'
    'model { '
    '  //priors'
    '  rho ~ uniform(-1,1);'
    '  mu ~ normal(0,inv_sqrt(0.001));    //uninformative prior on mean'
    '  tau ~ gamma(0.001,0.001);          //uninformative prior on precision'
    '  //likelihood'
    '  x ~ multi_normal(mu,E);            //MVN likelihood'
    '}'
    'generated quantities {'
    '  real prior_rho;                    //prior predictive'
    '  prior_rho = uniform_rng(-1,1);'
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
fit =  sm.sampling('file',  stanFilePath, ...
                   'model_name',   modelName, ...
                   'sample_file',  modelName, ...
                   'verbose',      true, ...
                   'chains',       nChains, ...
                   'warmup',       nWarmup, ...
                   'iter',         nIterations, ...
                   'working_dir',  workingDir);
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
tracedensity(samples,{'mu','rho','sigma'},diagnostics)
rankplots(samples,{'mu','rho','sigma'})

%% other plots
%plot densities
multidensity(samples,{'rho','mu','sigma'})

%% parameter recovery
estimatedValues = getsamplestats(samples,trueValues);
recoveryCounts = ...
    plotrecovery(estimatedValues,trueValues,{'sigma','rho','mu'});
proportionRecovered = recoveryCounts(1)/sum(recoveryCounts);
fprintf('\n%i of %i parameters (%.2f%%) were sucessfully recovered!\n', ...
    recoveryCounts(1),sum(recoveryCounts),proportionRecovered*100)

%% hypothesis test
plotdensity(samples,'rho','mean','credible95','credible90','credible50', ...
    trueValues.rho,'legend')
xlim([-1 1])

fprintf('\nthe true pearson correlation (rho) was: %.3f\n', ...
    trueValues.rho)
fprintf('the point estimate of the pearson correlation is: %.3f\n', ...
    estimatedValues.rho.mean)
if estimatedValues.rho.inCI95; inCI = ''; else; inCI = ' NOT'; end
fprintf(['this estimate for rho is%s in the 95%% credible interval: ' ...
    '[%.3f %.3f]\n'],inCI,estimatedValues.rho.lowerCI95, ...
    estimatedValues.rho.upperCI95)

%% hypothesis test
compareAt = 0;
[BF10,BF01] = savagedickey( ...
    samples.prior_rho,samples.rho,compareAt,'support',[-1 1]);
if BF10 > 1
    fprintf(['\nthe savage dickey bayes factor supports the ' ...
        'alternative: BF10 = %.3g\n'],BF10)
else
    fprintf(['\nthe savage dickey bayes factor supports the ' ...
        'null: BF01 = %.3g\n'],BF01)
end
