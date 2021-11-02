clear
clc

%in this model, a participant completes an N-arm bandit problem (i.e., 
%chooses among N response options).  rewards are stochastic.
%individual differences are accomodated by including hyperpriors in the 
%model. 

project = 'RL_bandit';

%% inputs

%experiment settings
RewardProb = [0.8 0.2];     %actual probability of reward per arm
                            %(length of vector == number of arms)
                            %(elements must sum to 1)

%mcmc settings
nChains     = 4;
nBurn       = 10;
nIterations = 1000;
nThin       = 1;

%an absolute path to a location for saved output
% output_dir = '/my/custom/dir/';
% output_dir = pwd;
outputDir = [pwd filesep 'output_' project];

%an absolute path to a location for Stan output files
% working_dir = '/my/custom/dir/';
% working_dir = output_dir;
workingDir = ['/tmp/wdir_' project];

%delete all temp files after stan finishes and selected output is saved?
clean_up_working_dir = true;

%% prepare for output
fprintf('\n\n**********\n\npreparing to run the %s model.\n', ...
    upper(project))

%output directory
if ~isequal(outputDir(end),filesep)
    outputDir(end+1) = filesep;
end
fprintf('\noutput will be saved to:\n> %s\n',outputDir)
if ~exist(outputDir,'dir')
    fprintf(['    currently this directory does not exist.  \n' ...
        '    making the directory now ... '])
    mkdir(outputDir)
    fprintf('done.\n')
end

%working directory
if ~isequal(workingDir(end),filesep)
    workingDir(end+1) = filesep;
end
fprintf(['\nStan files will be saved to the working directory:\n' ...
    '> %s\n'], workingDir)
if ~isequal(outputDir,workingDir)
    if ~exist(workingDir,'dir')
        fprintf(['    currently this directory does not exist.  \n' ...
            '    making the directory now ... '])
        mkdir(workingDir)
        fprintf('done.\n')
    end
end

%% data
%%% true parameter values
fprintf('\nsimulating data ... ')
%known features of the experiment
nSubjects = 20;                             %number of subjects
nTrials = randi([20 22],[nSubjects 1]);     %number of trials per subject
nData = sum(nTrials);                       %number of data points

nArms = length(RewardProb);                 %number of arms in bandit problem
bestArm = find(RewardProb == max(RewardProb)); %highest probability arm

%fixed parameters
Q0 = ones([nArms 1])/nArms;
epsilon = 0.001;

%enforced parameter value limits
%***included as distribution truncations in model specification***
b_limits = [0 40];
sig_b_limits = [0 20];
a_limits = [0 1];
sig_a_limits = [0 0.5];

%true parameter values --- group level
%beta
mu_beta = 8;
sigma_beta = 2;
%alpha  %DO NOT CHANGE FORM OF PRIOR --- YOU WILL HAVE PROBLEMS
a_alpha = 1.5;
b_alpha = 10;
mu_alpha = a_alpha/(a_alpha + b_alpha);
sigma_alpha = sqrt(a_alpha*b_alpha/((a_alpha+b_alpha)^2)/(a_alpha+b_alpha+1));

%true parameter values --- participant level
beta = abs(normrnd(mu_beta,sigma_beta,[nSubjects 1]));
alpha = betarnd(a_alpha,b_alpha,[nSubjects 1]);

%create a structure of true values
trueValues = collecttruevalues(beta,alpha, ...
    mu_beta,sigma_beta,mu_alpha,sigma_alpha);
disp(trueValues)

%%% simulated observations
%data IDs
Subject = NaN([nData 1]);
Trial = NaN([nData 1]);
%data
Action = NaN([nData 1]);
Reward = NaN([nData 1]);
Correct = NaN([nData 1]);

softmax = @(x) exp(x)/sum(exp(x));

n = 0; %(initialize data point counter)
for s = 1:nSubjects
    for t = 1:nTrials(s)
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
        Trial(n) = t;
        Action(n) = action;
        Reward(n) = r;
        Correct(n) = action == bestArm;
    end
end

%%%%%
% whos nSubjects nTrials nData Subject Trial Action Reward
% disp([Subject,Trial,Action,Reward])

%collect data for Stan
data = struct('NS',nSubjects,'NT',nTrials','NX',nData,'NA',nArms, ...
    'Subject',Subject,'Trial',Trial,'Action',Action,'Reward',Reward, ...
    'Q0',Q0);

fprintf('done!\n')

%% plot data
%group-level learning curve
figure('color',[1 1 1])
hold on
plot([0 min(nTrials)+1],(1/nArms)*[1 1],'--','color',getcolors('gray'));
accMean = NaN([min(nTrials) 1]);
accSE = NaN([min(nTrials) 1]);
for t = 1:min(nTrials)
    nthTrialCorr = Correct(Trial == t);
    accMean(t) = mean(nthTrialCorr);
    accSE(t) = std(nthTrialCorr)/sqrt(length(nthTrialCorr));
end
errorbar(accMean,accSE,'k')
%format plot
ylim([0 1])
xlabel('trial')
ylabel('accuracy')
set(gca,'box','on')

%choice density plot
f = figure('color',[1 1 1]);
fpos = get(f,'pos');
set(f,'pos',[fpos(1:2) 250 + 10*min(nTrials) 100+40*nArms])
hold on
choiceProportion = NaN([min(nTrials) nArms]);
choiceSE = NaN([min(nTrials) 1]);
for t = 1:min(nTrials)
    for a = 1:nArms
        choiceProportion(t,a) = sum(Action(Trial == t) == a)/nSubjects;
    end
end
choiceColors = getcolors('pastel');
for t = 1:min(nTrials)
    [~,preferredArmThisTrial] = max(choiceProportion(t,:));
    for a = 1:nArms
        h = plot(t,a,'linestyle','none', ...
        'marker','s','markersize',choiceProportion(t,a)*10+0.01, ...
        'markeredgecolor',choiceColors(a,:),'markerfacecolor','none');
        if a == preferredArmThisTrial
            h.MarkerFaceColor = choiceColors(a,:);
        end
    end
end
%format plot
set(gca,'ylim',[0.5 nArms+0.5],'ytick',1:nArms)
xlabel('trial')
ylabel('choice')

pause(1) %to let figures load

%% model specification
sp = @sprintf;

model_code = {
    'data { '
    '    int<lower=1> NS;                       //number of subjects'
    '    int<lower=1> NT[NS];                   //number of trials per subject'
    '    int<lower=1> NX;                       //number of data points'
    '    int<lower=1> NA;                       //number of bandit arms'
    '    int<lower=1,upper=NS> Subject[NX];     //subject number'
    '    int<lower=1,upper=max(NT)> Trial[NX];  //trial number'
    '    int<lower=1,upper=NA>  Action[NX];     //action selected'
    '    int<lower=0,upper=1>  Reward[NX];      //reinforcement given'
    '    vector<lower=0,upper=1>[NA] Q0;        //inital value of actions'
    '}'
    'parameters { '
    '    //group-level parameters'
    '    real<lower=0> shape_beta;'
    '    real<lower=0> scale_beta;'
%  sp('    real<lower=%.1f,upper=%.1f> mu_beta;',b_limits)
%  sp('    real<lower=%.1f,upper=%.1f> sigma_beta;',sig_b_limits)
 sp('    real<lower=%.2f,upper=%.2f> mu_alpha;',a_limits)
 sp('    real<lower=%.2f,upper=%.2f> sigma_alpha;',sig_a_limits)
    '    '
    '    //subject-level parameters'
 sp('    real<lower=%.1f,upper=%.1f> beta[NS];',b_limits)
 sp('    real<lower=%.2f,upper=%.2f> alpha[NS];',a_limits)
    '}'
    'transformed parameters {'
    '  real mu_beta = shape_beta/scale_beta;'
    '  real sigma_beta = sqrt(shape_beta / (scale_beta^2));'
    '}'
    'model { '
    '    //local variables'
    '    vector[NA] Q;'
    '    vector[NA] pi;'
    '    real delta;'
    '    real epsilon = 0.001;'
    '    '
    '    //group-level priors'
    '    shape_beta ~ gamma(2,2);'
    '    scale_beta ~ gamma(2,2);'
%  sp('    mu_beta ~ uniform(%.1f,%.1f);',b_limits)
%  sp('    sigma_beta ~ gamma(2,0.5) T[%.1f,%.1f];',sig_b_limits)
 sp('    mu_alpha ~ uniform(%.2f,%.2f);',a_limits)
 sp('    sigma_alpha ~ uniform(%.2f,%.2f);',sig_a_limits)
    '    '
    '    //subject-level priors'
    '    for (s in 1:NS) {'
%  sp('    beta[s] ~ normal(mu_beta,sigma_beta) T[%.1f,%.1f];',b_limits)
    '    beta[s] ~ gamma(shape_beta,scale_beta);'
 sp('    alpha[s] ~ normal(mu_alpha,sigma_alpha) T[%.2f,%.2f];',a_limits)
    '    }'
    '    '
    '    //likelihood'
    '    for (n in 1:NX) { '
    '        if (Trial[n]==1) {'
    '          for (a in 1:NA) {'
    '            Q[a] = Q0[a];'
    '          }'
    '        }'
    '        '
    '        //choice behavior'
    '        pi = (1-epsilon)*softmax(beta[Subject[n]]*Q) + epsilon/NA;'
    '        Action[n] ~ categorical(pi);'
    '        '
    '        //reward'
    '        delta = Reward[n]-Q[Action[n]];'
    '        Q[Action[n]] <- Q[Action[n]] + alpha[Subject[n]]*delta;'
    '    }'
    '    '
    '}'
};

%write the model code to a .stan file
stanfile = [workingDir project '.stan'];
stanfileID = fopen(stanfile,'w');
fprintf(stanfileID,'%s\n',model_code{:});
fclose(stanfileID);
%copy file to output directory
copyfile(stanfile,outputDir)

%% compile the model
tic
fprintf('\ncompiling the model ... ')
sm = StanModel('file',stanfile);
sm.compile();
fprintf('done!\n')
runtime = toc;
fprintf('compiling took %.2f seconds.\n',runtime)
    
%% run the model
tic
fprintf('\nrunning the model ... \n\n**********\n\n')
% sm.control.adapt_delta = 0.9;
fit =  sm.sampling('file',  stanfile, ...
            'model_name',   project, ...
            'sample_file',  project, ...
            'data',         data, ...
            'chains',       nChains, ...
            'warmup',       nBurn, ...
            'iter',         nIterations, ...
            'thin',         nThin, ...
            'verbose',      true, ...
            'working_dir',  workingDir);
fit.block();
stanSummary = fit.print;
fprintf('\n**********\n\ndone!\n')

runtime = toc;
fprintf('\nsampling took %.2f seconds.\n',runtime)

%% extract from the StanFit object
[samples,diagnostics] = extractsamples('matlabstan',fit);
parameters = fieldnames(samples);
instances = getparaminstances([],samples);

%clean up after MATLABStan
clearvars fit
%clean up after Stan
if clean_up_working_dir, delete([workingDir '*']); rmdir(workingDir); end

%% diagnostic report
%calculate any convergence diagnostics based on posterior samples
posteriorTable = mcmctable(samples);
%print a report about all mcmc diagnostics
interpretdiagnostics(diagnostics,posteriorTable)

%% get parameter estimates & other statistics
estimatedValues = getsamplestats(samples,trueValues);

%% output reports & plots
figure;
plotlp(diagnostics,posteriorTable)

recoveryCounts = [0 0];
requests = {'alpha','beta', ...
    {'mu_alpha','sigma_alpha','mu_beta','sigma_beta'}};
for n = 1:length(requests)
    rc = plotrecovery(estimatedValues,trueValues,requests{n});
    recoveryCounts = recoveryCounts + rc;
end
proportionRecovered = recoveryCounts(1)/sum(recoveryCounts);
fprintf('\n%i of %i parameters (%.2f%%) were sucessfully recovered!\n', ...
    recoveryCounts(1),sum(recoveryCounts),proportionRecovered*100)

%% save figures
%save figures? where? 
figFolder = [pwd filesep 'figures'];
figFormats = {'pdf','png'};%{'pdf','png'};

%figure saving
if figFolder(end)~=filesep, figFolder(end+1) = filesep; end
if ~isfolder(figFolder), mkdir(figFolder); end
for ff = 1:length(figFormats)
    if figFormats{ff}(1)~= '.', figFormats{ff} = ['.' figFormats{ff}]; end
%     figFormats{ff} = ['_' figFormats{ff}];
end

figure % (ensures MATLAB doesn't goof the sizing on the first plot)
title('dummy figure (to protect figure sizing)')

%%%

jointdensity(samples,'sigma_alpha','mu_alpha',diagnostics)
figName = 'RL___jointdensity';
for ff = 1:length(figFormats), exportgraphics(gcf,[figFolder figName figFormats{ff}]); end

multidensity(samples,{'mu_alpha','sigma_alpha','mu_beta','sigma_beta'},diagnostics)
figName = 'RL___multidensity';
for ff = 1:length(figFormats), exportgraphics(gcf,[figFolder figName figFormats{ff}]); end

tracedivergent(samples,'sigma_alpha',diagnostics)
figName = 'RL___tracediv';
for ff = 1:length(figFormats), exportgraphics(gcf,[figFolder figName figFormats{ff}]); end

parcoordivergent(samples,diagnostics,{'mu_alpha','sigma_alpha','alpha[1]','alpha[2]','alpha[3]'})
figName = 'RL___parcoor';
for ff = 1:length(figFormats), exportgraphics(gcf,[figFolder figName figFormats{ff}]); end

plotlp(diagnostics,rtable)
figName = 'RL___plotlp';
for ff = 1:length(figFormats), exportgraphics(gcf,[figFolder figName figFormats{ff}]); end
