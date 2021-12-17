
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%      EXAMPLE_RL       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script implements a hierarchical cogntive model --- specifically a 
% simple reinforcement learning model of data from an N-armed bandit task.  
% 
% first, data from this experimental design is simulated for multiple 
% participants using true parameter values (randomly drawn from the priors).  
% then, a Bayesian implementation of Q-learning is specified in Stan code.
% finally, the model is run and parameter recovery is demonstrated. 
% 
% this example is designed to illustrate the following ideas & principles:
%   - implementation of a hierarchical cognitive model
%   - computational tricks for beta and gamma priors & their associated
%       hyperpriors
%   - ...
% 
% 
% to run this script, the following must be installed: 
%       Stan        >>  https://mc-stan.org/
%       MATLABStan  >>  https://github.com/brian-lau/MatlabStan
%       matstanlib  >>  https://github.com/baribault/matstanlib
% 
% Reference:    ... 
% 
% (c) beth baribault 2021 ---                                 > matstanlib
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% inputs

%which model specification to use (for both simulation & fitting)?
% modelName = 'RL_broken'; %uses normal & uniform priors
modelName = 'RL_fixed';  %uses beta & gamma priors

%experimental design
%what is the actual probability of reward for each bandit arm?
%   (specify this as a vector of reward probabilites where the 
%    number of elements == number of arms, and the elements sum to 1)
% RewardProb = [1 0 0]; %deterministic, 3-armed bandit
RewardProb = [0.7 0.1 0.2]; %probabilistic, 3-armed bandit

%sampler settings
nChains     = 4;        %how many chains?
nWarmup     = 250;      %how many iterations during warmup?
nIterations = 100;     %how many iterations after warmup?

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

%double check RewardProb
if ~isnumeric(RewardProb) || ~isvector(RewardProb) || length(RewardProb) < 2
    error('RewardProb be a numeric vector with at least two elements.')
elseif ~isequal(sum(RewardProb),1)
    error('RewardProb must sum to 1.')
end

%% simulate data
fprintf('\nsimulating data ... \n')

%%% known features of the experiment
nSubjects = 30;     %number of subjects
nProblems = 4;      %number of bandit problems per subject
nTrials   = 20;     %number of trials per problem

nArms = length(RewardProb); %number of arms in bandit problem
bestArm = find(RewardProb == max(RewardProb)); %highest probability arm
nData = nSubjects*nProblems*nTrials; %total number of data points


%%% fixed parameters
% Q0 = ones([nArms 1])/nArms;
Q0 = zeros([nArms 1]);
epsilon = 0.005;

%%% true parameter values
switch modelName
    case 'RL_fixed'
        %beta
        mu_beta = 7.5; b1plus = 5; 
        b1 = b1plus-1; b2 = b1plus/mu_beta;
        sigma_beta = sqrt(b1plus/b2^2);
        beta = gamrnd(b1plus,1/b2,[nSubjects 1]);
        %alpha
        a1plus = 3; a2plus = 6;         %actually used in the prior
        a1 = a1plus-1; a2 = a2plus-1;   %inferred 
        mu_alpha = a1plus/(a1plus+a2plus);
        sigma_alpha = sqrt(a1plus*a2plus / ((a1plus+a2plus)^2 * (a1plus+a2plus+1)) );
        alpha = betarnd(1+a1,1+a2,[nSubjects 1]);
        %phi
        p1plus = 1.5; p2plus = 10;      %actually used in the prior
        p1 = p1plus-1; p2 = p2plus-1;   %inferred
        mu_phi = p1plus/(p1plus+p2plus);
        sigma_phi = sqrt(p1plus*p2plus / ((p1plus+p2plus)^2 * (p1plus+p2plus+1)) );
        phi = betarnd(1+p1,1+p2,[nSubjects 1]);
    case 'RL_broken'
        %beta
        mu_beta = 7.5; sigma_beta = 2.5;
        bCheck = false; 
        while ~bCheck
            beta = normrnd(mu_beta,sigma_beta,[nSubjects 1]);
            if all(beta > 0), bCheck = true; end
        end
        %alpha
        mu_alpha = 0.25; sigma_alpha = 0.15;
        aCheck = false; 
        while ~aCheck
            alpha = normrnd(mu_alpha,sigma_alpha,[nSubjects 1]);
            if all(alpha > 0 & alpha < 1), aCheck = true; end
        end
        %phi
        mu_phi = 0.075; sigma_phi = 0.025;
        pCheck = false; 
        while ~pCheck
            phi = normrnd(mu_phi,sigma_phi,[nSubjects 1]);
            if all(phi > 0 & phi < 1), pCheck = true; end
        end
    otherwise
        error('''%s'' is not a recognized modelName input.',modelName)
end

%%% simulate behavior
%data IDs
Subject = NaN([nData 1]);
Problem = NaN([nData 1]);
Trial = NaN([nData 1]);
%data
Action = NaN([nData 1]);
Reward = NaN([nData 1]);
Correct = NaN([nData 1]);

%introduce softmax as an anonymous function 
softmax = @(x) exp(x)/sum(exp(x));

n = 0; %initialize data point counter
for s = 1:nSubjects
  for p = 1:nProblems
    for t = 1:nTrials
        %if first trial, initialize value
        if t == 1
            Q = Q0;
        else
            Q = (1-phi(s))*Q + phi(s)*Q0;
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
trueValues = collecttruevalues(beta,alpha,phi, ...
    mu_beta,sigma_beta,mu_alpha,sigma_alpha,mu_phi,sigma_phi);
disp(trueValues)

%collect data for Stan
dataStruct = struct('S',nSubjects,'T',nTrials','A',nArms,'N',nData, ...
    'Subject',Subject,'Trial',Trial,'Action',Action,'Reward',Reward, ...
    'Q0',Q0);

fprintf('done!\n')

%% plot the simulated behavior
fprintf('plotting the simulated behavioral data ...')

msl.options %load plotting options

dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) [420 280]*figScaling];
close(dumf.Number); %close dummy figure
hold on

%%%%%%%%%% data
%subject-level learning curves
nTrialsPerBin = 2;
nTrialBins = floor(nTrials/nTrialsPerBin);
trialBinCenters = (1:nTrialBins)*nTrialsPerBin - 0.5*(nTrialsPerBin-1);

subColors = makecolormap(0.775*[1 1 1],0.95*[1 1 1],nSubjects);
subMeans = computeLearningCurve( ...
    table(Trial,Correct,Subject),'Correct','Trial', ...
    {'Subject'},[],nTrialsPerBin);
for s = 1:nSubjects
    plot(trialBinCenters,subMeans(:,s)','color',subColors(s,:),'linewidth',1)
end

%group-level learning curve
nTrialsPerBin = 1;
nTrialBins = floor(nTrials/nTrialsPerBin);
binLabels = ...
    cellfun(@num2str,num2cell(nTrialsPerBin:nTrialsPerBin:nTrials),'uni',0);
trialBinCenters = (1:nTrialBins)*nTrialsPerBin - 0.5*(nTrialsPerBin-1);

[groupMean,groupSEM] = computeLearningCurve( ...
    table(Trial,Correct,Subject),'Correct','Trial', ...
    {},[],nTrialsPerBin);
errorbar(trialBinCenters,groupMean,groupSEM,'k')
plot(trialBinCenters,groupMean,'color','k','linewidth',2.5)

%overlay line at chance
plot([0 min(nTrials)+1],(1/nArms)*[1 1],':','color','k');

%format plot
xlim([0.5 nTrials+0.5])
set(gca,'xtick',4:4:nTrials)
ylim([0 1.005])
xlabel('trial')
ylabel('proportion correct')
set(gca,'box','on','fontsize',fontSz)

pause(0.1) %to let figures load
fprintf('done!\n')

%% compile the model

%copy the .stan file containing the model code to the working directory
stanFile = [modelName '.stan'];
stanFilePath = [workingDir stanFile];
copyfile(which([modelName '.stan']),workingDir)

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

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if cleanUp, delete([workingDir '*']); rmdir(workingDir); end

%% diagnostic reports & plots
%extract parameter names
parameters = {'beta','alpha','phi'};
hyperparameters = getparaminstances('*_*',samples);

%compute posterior sample-based diagnostics and summary statistics
posteriorTable = mcmctable(samples);
%print a report about all MCMC diagnostics
interpretdiagnostics(diagnostics,posteriorTable)

% %trace plots/rank plots
% tracedensity(samples,hyperparameters,diagnostics)
% rankplots(samples,hyperparameters)

multidensity(samples,{'*_b*','*_p*'},diagnostics,{'divergent','treedepth','energy'})

parallelsamples(samples,diagnostics, ...
    {'mu_alpha','sigma_alpha','mu_phi','sigma_phi','alpha[1]','phi[1]'})
parallelsamples(samples,diagnostics, ...
    {'mu_beta','sigma_beta','mu_phi','sigma_phi','mu_alpha','sigma_alpha'},true)

%% get parameter estimates & other statistics

%plot recovery
recoveryCounts = [0 0];
recoveryCounts = recoveryCounts + ...
    plotrecovery(samples,trueValues,'beta', ...
    'addinterval',true);
recoveryCounts = recoveryCounts + ...
    plotrecovery(samples,trueValues,'alpha', ...
    'addinterval',true);
recoveryCounts = recoveryCounts + ...
    plotrecovery(samples,trueValues,'phi', ...
    'addinterval',true);
recoveryCounts = recoveryCounts + ...
    plotrecovery(samples,trueValues,hyperparameters);
%report recovery
proportionRecovered = recoveryCounts(1)/sum(recoveryCounts);
fprintf('\n%i of %i parameters (%.2f%%) were sucessfully recovered!\n', ...
    recoveryCounts(1),sum(recoveryCounts),proportionRecovered*100)

%% visualize group-level means
% % plotdensity(samples,'mu_alpha','credible','mean','legend',mean(trueValues.alpha))
% % xlim([0 1])
% % ylabel('')
% % 
% % plotdensity(samples,'mu_beta','credible','hdi','median','legend',mean(trueValues.beta))
% % xlim([0 17])
% % ylabel('')

%% visualize individual differences
% plotintervals(samples,'beta', ...
% 	'truevalues',trueValues,'truecolor',getcolors('blue'))
% xlim([0 20])
% % % horzdensity(samples,'beta', ...
% % %     'truevalues',trueValues,'truecolor',getcolors('blue'))
% 
% plotintervals(samples,'alpha', ...
%     'truevalues',trueValues,'truecolor',getcolors('blue'))
% xlim([0 1])
% % % horzdensity(samples,'alpha', ...
% % %     'truevalues',trueValues,'truecolor',getcolors('blue'))
% 
% plotintervals(samples,'phi', ...
% 	'truevalues',trueValues,'truecolor',getcolors('blue'))
% xlim([0 1])
% % % horzdensity(samples,'phi', ...
% % %     'truevalues',trueValues,'truecolor',getcolors('blue'))
