
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script demonstates how a hierarchical cognitive model may be
% implemented in Stan.  
% 
% specifically, it demonstrates a reinforcement learning model of data from
% an N-armed bandit task. first, data from this experimental design is
% simulated for multiple participants.  then, a Bayesian implementation of
% Q-learning is specified as a Stan model.  finally, the model is run and
% parameter recovery is demonstrated. 
% 
% this example is designed to illustrate the following ideas & principles:
%   - implementation of a hierarchical cognitive model
%   - computational tricks for beta and gamma priors & their associated
%       hyperpriors
%   - 
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib 
%
% 
% (c) beth baribault 2021 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

modelName = 'example_RL';

%experiment settings
RewardProb = [0.8 0.2];     %actual probability of reward per arm
                            %(length of vector == number of arms)
                            %(elements must sum to 1)

%which model specification, 'gamma' priors or truncated 'normal' priors?
priors = 'gamma';

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 50;%100;      %how many iterations during warmup?
nIterations = 200;%500;      %how many iterations after warmup?

%an absolute path to a location for (temporary) Stan output files
% workingDir = '/my/custom/dir/';
workingDir = [pwd filesep 'wdir_' modelName];

%delete all temp files after stan finishes and selected output is saved?
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

%%% known features of the experiment
nSubjects = 25;     %number of subjects
nProblems = 4;      %number of bandit problems per subject
nTrials   = 12;     %number of trials per bandit problem
nData = nSubjects*nProblems*nTrials; %total number of data points

nArms = length(RewardProb); %number of arms in bandit problem
bestArm = find(RewardProb == max(RewardProb)); %highest probability arm

%%% fixed parameters
Q0 = ones([nArms 1])/nArms;
epsilon = 0.001;

%%% true parameter values
% %... participant level
% beta = gamrnd(shape,scale,[nSubjects 1]);
% alpha = betarnd(a,b,[nSubjects 1]);

% beta = linspace(1,20,nSubjects);
% alpha = ones([1 nSubjects])*0.25;

%%% true parameter values
N1 = floor(nSubjects/2);
N2 = nSubjects - N1;

%... participant level
beta  = [rand([1 N1])*5 + 5         linspace(5,10,N2)         ];
alpha = [linspace(0.05,0.25,N1)     rand([1 N2])*0.225 + 0.005];

%%% simulate observations
%data IDs
Subject = NaN([nData 1]);
Problem = NaN([nData 1]);
Trial = NaN([nData 1]);
%data
Action = NaN([nData 1]);
Reward = NaN([nData 1]);
Correct = NaN([nData 1]);

softmax = @(x) exp(x)/sum(exp(x));

n = 0;
for s = 1:nSubjects
  for p = 1:nProblems
    for t = 1:nTrials
        %if first trial, initialize value
        if t == 1
            Q = Q0;
        end
        %simulate action
        p_action = (1-epsilon)*softmax(beta(s)*Q) + epsilon/nArms;
        action = find(mnrnd(1,p_action));
        %simulate reinforcement
        r = binornd(1,RewardProb(action));
        Q(action) = Q(action) + alpha(s)*(r - Q(action));
        
        %record simulated data
        n = n + 1; %(increment data point counter)
        Subject(n) = s;
        Problem(n) = p;
        Trial(n) = t;
        Action(n) = action;
        Reward(n) = r;
        Correct(n) = action == bestArm;
    end
  end
end

%%%%%
% whos nSubjects nTrials nData Subject Trial Action Reward
% disp([Subject,Trial,Action,Reward])

%create a structure of true values for observing parameter recovery
trueValues = collecttruevalues(beta,alpha);
disp(trueValues)

%collect data for Stan
data = struct('S',nSubjects,'T',nTrials','N',nData,'A',nArms, ...
    'Subject',Subject,'Trial',Trial,'Action',Action,'Reward',Reward, ...
    'Q0',Q0);

disp('done!')

%% plot data
matstanlib_options %load plotting options

dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) [420 280]*figScaling];
close(dumf.Number); %close dummy figure
hold on

%%%%%%%%%% data
%subject-level learning curves
nTrialsPerBin = 2;
nBins = floor(nTrials/nTrialsPerBin);
binLocs = (1:nBins)*nTrialsPerBin - (nTrialsPerBin-1)/2;

subColors = makecolormap(0.775*[1 1 1],0.95*[1 1 1],nSubjects);
subMeans = computeLearningCurve( ...
    table(Trial,Correct,Subject),'Correct','Trial', ...
    {'Subject'},[],nTrialsPerBin);
for s = 1:nSubjects
    plot(binLocs,subMeans(:,s)','color',subColors(s,:),'linewidth',1)
end

%group-level learning curve
nTrialsPerBin = 1;
nBins = floor(nTrials/nTrialsPerBin);
binLabels = cellfun(@num2str,num2cell(nTrialsPerBin:nTrialsPerBin:nTrials),'uni',0);
binLocs = (1:nBins)*nTrialsPerBin - (nTrialsPerBin-1)/2;

[groupMean,groupSEM] = computeLearningCurve( ...
    table(Trial,Correct,Subject),'Correct','Trial', ...
    {},[],nTrialsPerBin);
errorbar(binLocs,groupMean,groupSEM,'k')
plot(binLocs,groupMean,'color','k','linewidth',2.5)

%overlay line at chance
plot([0 min(nTrials)+1],(1/nArms)*[1 1],':','color','k');

%format plot
xlim([0.5 nTrials+0.5])
set(gca,'xtick',4:4:nTrials)
ylim([0.25 1.01])
xlabel('trial')
ylabel('proportion correct')
set(gca,'box','on','fontsize',fontSz)

pause(1) %to let figures load
clearvars f

%% model specification

switch priors
  case 'gamma'
    %using a gamma prior for beta and a beta prior for alpha, and wide
    %gamma priors for all hyperparameters
    hyperpriortester([0 25;0 15], ...
        'gamma', ...            %prior distribution
        {5,1.5}, ...            %ideal prior
        'gamma',{10,1/2}, ...   %hyperprior for shape
        'gamma',{3,1/2}, ...    %hyperprior for rate
        {@(x) x+1, @(x) x})     %functions applied to each hyperparameter
    hyperpriortester([0 1;0 15], ...
        'beta', ...            %prior distribution
        {1.1,1.1}, ...         %ideal prior
        'gamma',{1,1/1}, ...   %hyperprior for shape
        'gamma',{1,1/1}, ...   %hyperprior for rate
        {@(x) x+1, @(x) x+1})  %functions applied to each hyperparameter
    pause(0.1)
    
    modelCode = {
        'data { '
        '  int<lower=1> S;                      //number of subjects'
        '  int<lower=1> T;                      //number of trials per subject'
        '  int<lower=1> N;                      //number of data points'
        '  int<lower=1> A;                      //number of bandit arms'
        '  int<lower=1,upper=S> Subject[N];     //subject number'
        '  int<lower=1,upper=T> Trial[N];       //trial number'
        '  int<lower=1,upper=A> Action[N];      //action selected'
        '  int<lower=0,upper=1> Reward[N];      //reinforcement given'
        '  vector[A] Q0;                        //inital value of actions'
        '}'
        'transformed data { '
        '  //fixed parameters'
        '  real epsilon = 0.001;'
        '}'
        'parameters { '
        '  //group-level parameters'
        '  real<lower=0> shape;'
        '  real<lower=0> rate;'
        '  real<lower=0> a;'
        '  real<lower=0> b;'
        '  '
        '  //subject-level parameters'
        '  real<lower=0> beta[S];               //inverse temperature'
        '  real<lower=0,upper=1> alpha[S];      //learning rate'
        '}'
        'transformed parameters { '
        '  //derived parameters'
%         '  real mu_beta = (shape+1)/rate;'
%         '  real sigma_beta = sqrt((shape+1) / (rate^2));'
%         '  real mu_alpha = (a+1)/(a+1 + b+1);'
%         '  real sigma_alpha = sqrt((a+1)*(b+1)/( ((a+1 + b+1)^2)*(a+1 + b+1 + 1) ));'
        '  real mu_beta = shape/rate;'
        '  real sigma_beta = sqrt(shape / (rate^2));'
        '  real mu_alpha = a/(a + b);'
        '  real sigma_alpha = sqrt(a*b/( ((a+b)^2)*(a+b+1) ));'
        '}'
        'model { '
        '  //local variables'
        '  vector[A] Q;'
        '  vector[A] pi;'
        '  real delta;'
        '  '
        '  //group-level priors'
        '  shape ~ gamma(1,1);'%(10,3);'
        '  rate ~ gamma(1,1);';%(3,2);'
        '  a ~ gamma(1,1);'
        '  b ~ gamma(1,1);'
        '  '
        '  //individual-level priors'
        '  for (s in 1:S) { '
        '    beta[s] ~ gamma(shape,rate);'%1+shape,rate);'
        '    alpha[s] ~ beta(a,b);'%1+a,1+b);'
        '  }'
        '  '
        '  //likelihood'
        '  for (n in 1:N) { '
        '    //initialize value'
        '    if (Trial[n]==1) { '
        '      Q = Q0;'
        '    }'
        '    '
        '    //choice behavior'
        '    pi = (1-epsilon)*softmax(beta[Subject[n]]*Q) + epsilon/A;'
        '    Action[n] ~ categorical(pi);'
        '    '
        '    //reward'
        '    delta = Reward[n] - Q[Action[n]];'
        '    Q[Action[n]] = Q[Action[n]] + alpha[Subject[n]]*delta;'
        '  }'
        '  '
        '}'
    };
  case 'normal'
    %using normal priors on alpha/beta, and uniform priors for
    %the hyperparameters
    modelCode = {
        'data { '
        '  int<lower=1> S;                      //number of subjects'
        '  int<lower=1> T;                      //number of trials per subject'
        '  int<lower=1> N;                      //number of data points'
        '  int<lower=1> A;                      //number of bandit arms'
        '  int<lower=1,upper=S> Subject[N];     //subject number'
        '  int<lower=1,upper=T> Trial[N];       //trial number'
        '  int<lower=1,upper=A> Action[N];      //action selected'
        '  int<lower=0,upper=1> Reward[N];      //reinforcement given'
        '  vector[A] Q0;                        //inital value of actions'
        '}'
        'transformed data { '
        '  //fixed parameters'
        '  real epsilon = 0.001;'
        '}'
        'parameters { '
        '  //group-level parameters'
        '  real<lower=0,upper=1> mu_beta;'
        '  real<lower=0,upper=1> sigma_beta;'
        '  real<lower=0,upper=1> mu_alpha;'
        '  real<lower=0,upper=1> sigma_alpha;'
        '  '
        '  //subject-level parameters'
        '  real<lower=0> beta[S];               //inverse temperature'
        '  real<lower=0,upper=1> alpha[S];      //learning rate'
        '}'
        'model { '
        '  //local variables'
        '  vector[A] Q;'
        '  vector[A] pi;'
        '  real delta;'
        '  '
        '  //group-level priors'
        '  mu_alpha ~ uniform(0,40);'
        '  mu_beta ~ uniform(0,10);'
        '  mu_alpha ~ uniform(0,1);'
        '  sigma_alpha ~ uniform(0,0.5);'
        '  '
        '  //individual-level priors'
        '  for (s in 1:S) { '
        '    beta[s] ~ normal(mu_alpha,mu_beta);'
        '    alpha[s] ~ normal(mu_alpha,sigma_alpha);'
        '  }'
        '  '
        '  //likelihood'
        '  for (n in 1:N) { '
        '    //initialize value'
        '    if (Trial[n]==1) { '
        '      Q = Q0;'
        '    }'
        '    '
        '    //choice behavior'
        '    pi = (1-epsilon)*softmax(beta[Subject[n]]*Q) + epsilon/A;'
        '    Action[n] ~ categorical(pi);'
        '    '
        '    //reward'
        '    delta = Reward[n] - Q[Action[n]];'
        '    Q[Action[n]] = Q[Action[n]] + alpha[Subject[n]]*delta;'
        '  }'
        '  '
        '}'
    };
end

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
            'data',         data, ...
            'chains',       nChains, ...
            'warmup',       nWarmup, ...
            'iter',         nIterations, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
[stanSummaryTxt,stanSummary] = fit.print('sig_figs',5);
fprintf('\n**********\n\ndone!\n')
fprintf('\nsampling took %.2f seconds.\n',toc)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('MATLABStan',fit);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%% diagnostic report & plots
%calculate convergence diagnostics based on the posterior samples
% rtable = msl.rhattable(samples);
rtable = rhattable(samples);
%print a report about all mcmc diagnostics
interpretdiagnostics(diagnostics,rtable)

%extract parameter names
parameters = fieldnames(samples);
instances = getparaminstances([],samples);
hyperparameters = setdiff(parameters,{'alpha','beta'});

% %generate diagnostic plots
% % tracedensity(samples,hyperparameters,diagnostics)
% rankplots(samples,hyperparameters)
% 
% multidensity(samples,hyperparameters,diagnostics)
% 
% parcoordivergent(samples,diagnostics,{'mu_alpha','sigma_alpha', ...
%     'alpha[1]','alpha[2]','alpha[3]'})
% 
% if ismember('shape',hyperparameters)
%     jointdensity(samples,'shape','rate',diagnostics)
% else
%     jointdensity(samples,'mu_alpha','sigma_alpha',diagnostics)
% end

%% get parameter estimates & other statistics

estimatedValues = getsamplestats(samples,trueValues);

%plot recovery
recoveryCounts = [0 0];
recoveryCounts = recoveryCounts + plotrecovery(estimatedValues,trueValues,'alpha');
recoveryCounts = recoveryCounts + plotrecovery(estimatedValues,trueValues,'beta');
%report recovery
proportionRecovered = recoveryCounts(1)/sum(recoveryCounts);
fprintf('\n%i of %i parameters (%.2f%%) were sucessfully recovered!\n', ...
    recoveryCounts(1),sum(recoveryCounts),proportionRecovered*100)

%% visualize group-level means
plotdensity(samples,'mu_alpha','credible','mean','legend',mean(trueValues.alpha))
xlim([0 1])
ylabel('')

plotdensity(samples,'mu_beta','credible','hdi','median','legend',mean(trueValues.beta))
xlim([0 17])
ylabel('')

%% visualize individual differences
plotintervals(samples,'alpha','linecolor',getcolors('green'))
horzdensity(samples,'beta','fillcolor',getcolors('blue'))
