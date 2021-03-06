function hyperpriortester(xRange,priorForm,referenceParams,varargin)
% HYPERPRIORTESTER helps visualize hyperpriors' effects on a prior.
% 
% *** NOTE: MAKE SURE TO USE MATLAB'S DISTRIBUTION PARAMETERIZATIONS ***
% ***       WHEN CALLING THIS FUNCTION, THEN (if neccessary) CONVERT ***
% ***       TO STAN'S PARAMETERIZATIONS IN YOUR MODEL SPECIFICATION. ***
% 
% *** FOR SOME DISTRIBUTIONS, STAN & MATLAB USE DIFFERENT PARAMETERIZATIONS!
% 
% this function will test the effect of different hyperprior distributions
% on possible parameter values at the prior level, by visualizing how
% different hyperpriors distributions combine to constrain the prior. 
% 
% specifically, this function will generate a figure your candidate
% hyperprior distribution (with the given parameterization), simulate their
% joint distribution (to foresee challenges), and will contrast a reference
% prior with a simulated expected prior, based on the specified
% hyperpriors. 
% 
% NOTE: THIS FUNCTION IS HEAVILY RELIANT ON MATLAB'S STATS TOOLBOX.  IF YOU
% HAVE A QUESTION ABOUT THE ORDER IN WHICH PARAMETERS ARE EXPECTED, OR HOW
% A DISTRIBUTION IS PARAMETERIZED, PLEASE CONSULT THE MATLAB DOCUMENTATION.
% 
% *** REMEMBER: MATLAB AND STAN HAVE DIFFERENT PARAMETERIZATIONS ***
% ****          FOR SOME DISTRIBUTIONS!!! (e.g., Gamma)          ***
% 
% 
% HYPERPRIORTESTER(XRANGE,PRIORFORM,REFERENCEPARAMS, ...
%                  HYPERPRIORFORM,HYPERPRIORPARMS, ...
%                  [HYPERPRIORFORM2,HYPERPRIORPARMS2] [, ...])
%   XRANGE range of parameter values to consider (and plot!). 
%   if XRANGE is a two-element vector, then the same values will be used as
%   the domain to visualize both the prior and the hyperprior. 
%   alternatively, XRANGE may be a 2x2 matrix, where the first row is the
%   domain for the priors, and the second row is the domain for the
%   hyperpriors.  
%   
%   PRIORFORM is a string indicates what form of the prior you are using is (i.e., the
%   distribution name).   
%   
%   REFERENCEPARAMS contains parameter values which, in conjuction with
%   PRIORFORM, define a specific distribution that would be useful to
%   include as a reference.  for example, REFERENCEPARAMS might define an
%   ideal prior distribution (i.e., what prior one would use if the model
%   were not hierarchical.)  
%   if REFERENCEPARAMS is empty, no prior will be included as a
%   reference on the plot.
%   
%   for each parameter required by your prior's form, include one
%   HYPERPRIORFORM and one HYPERPRIORPARAM input.
%   HYPERPRIORFORM indicates what form of the [first] hyperprior you are
%   using is. 
%   
%   PRIORFORM and HYPERPRIORFORM should both be strings that are
%   able to be recognized by MATLAB's pdf function.  
%   a list (that may or may not be current given your MATLAB edition) is: 
%       'beta', 'binomial', 'chi-square', 'exponential', 'f', 
%       'gamma', 'gp' (generalized pareto), 'geometric', 'half normal', 
%       'hypergeometric', 'lognormal', 'nbin' (negative binomial), 
%       'noncentral f', 'noncentral t', 'noncentral chi-square', 
%       'normal', 'poisson', 'rayleigh', 't', 'unid' (discrete uniform), 
%       'unif' (continuous uniform), 'weibull'
%   type "help pdf" at the command line for a list of recognized strings &
%   distributions given YOUR edition of MATLAB, etc.
%   
%   PRIORPARAMS and HYPERPRIORPARAMS each must be cell type vectors where
%   all elements are scalar numeric values. 
% 
% 
% there are also a few optional inputs, which may be given in any order:
% 
% HYPERPRIORTESTER(...,'clf')
%   plots in the current figure (after wiping it), rather than creating a
%   new figure. 
% 
% HYPERPRIORTESTER(...,N)
%   bumps up the number of points in domain AND the number of simulations.
%   the default is 250.
% 
% HYPERPRIORTESTER(...,f)
%   f is a cell type vector of function handles of functions to apply to
%   each hyperparameter (that has been randomly drawn from the hyperpriors)
%   before the hyperparameter is used for prior simulation.   
%   (these functions can be used to, for example, shift and scale the
%   hyperpriors.  the only requirement is that the functions only have ONE
%   input: the parameter value!).
%   NOTE: these functions are not applied to the plotting of the
%   hyperpriors in the first panel.
% 
% 
% Example:
%   the following hyperpriors induce a distribtution over priors that
%   grossly overweights horseshoe-shaped priors (vs. the ideal prior). 
%   > hyperpriortester([0 1;0 10], ...       %prior & hyperprior domains
%                       'beta',{2,2}, ...    %an ideal prior
%                       'gamma',{1,1/1}, ... %hyperprior for A
%                       'gamma',{1,1/1})     %hyperprior for B
%   many of the priors place such extreme weight on 0 and 1 that adjusting
%   the y-axis limit is likely necessary to see what's going on: 
%   > ylim([0 8])
% 
%   setting a lower bound of 1 on both hyperparameters allows for a more
%   reasonable distribution of priors (that is, on average, closer to the
%   ideal prior, while still permitting a wide selection).  
%   > f = {@(x) x + 1; @(x) x + 1}; %add 1 to each hyperparameter before use
%   > hyperpriortester([0 1;0 10], ...       %prior & hyperprior domains
%                       'beta',{2,2}, ...    %an ideal prior
%                       'gamma',{1,1/1}, ... %hyperprior for A
%                       'gamma',{1,1/1}, ... %hyperprior for B
%                       f) %functions to be applied to the hyperparameters 
% 
% 
% See also PDF
% 
% (c) beth baribault 2020 ---                                 > matstanlib 

msl.options

%% parse inputs
nNamedInp = 3;
nVarargin = length(varargin);
if nargin < nNamedInp || nVarargin < 2
    error('too few inputs: at least four inputs total are required.')
end

%xrange
if isnumeric(xRange) && size(xRange,2)==2
    if any(xRange(:,1) >= xRange(:,2))
        error('xRange must be increasing (within each row).')
    end
else
    error(['the first input, xRange, must be a numeric ' ...
        'two-element vector (or two-column matrix).'])
end

%prior form
if ~ischar(priorForm)
    error('priorForm must be a string.')
end

%prior params
includeRef = true;
if isempty(referenceParams)
    includeRef = false;
elseif ~iscell(referenceParams) || ~isvector(referenceParams)
    error('REFERENCEPARAMS must be a cell type vector.')
elseif ~all(cellfun(@isnumeric,referenceParams) & cellfun(@isscalar,referenceParams))
    error('REFERENCEPARAMS must be a cell type vector of scalar numeric values.')
end

%prior form AGAIN
if ~ischar(priorForm)
    error('priorForm must be a string.')
else
    if includeRef
      try
        test = random(priorForm,referenceParams{:});
      catch
        error(['priorForm must be a valid distribution string.  ' ...
            'type "help pdf" for a list.'])
      end
      if isnan(test)
        error('invalid parameters given for prior form %s.',priorForm)
      end
    end
end

%% parse optional inputs

%initialize lists
hyperparamForms = {};
hyperparamParams = {};

%defaults
noNewFig = false;
NS = 250;

skipNext = false;
varargin{end+1} = NaN;
for v = 1:nVarargin
    if skipNext
        %skip parsing this input
        skipNext = false; %reset
    else
        %parse this input
        if isequal(varargin{v},'clf')
            noNewFig = true;
        elseif ischar(varargin{v}) && (iscell(varargin{v+1}) && isvector(varargin{v+1}))
            form = varargin{v};
            params = varargin{v+1};
            if all(cellfun(@isnumeric,params) & cellfun(@isscalar,params))
                try %see if valid distribution form name
                    test = random(form,params{:});
                catch 
                    error(['%s is not a valid distribution form name ' ...
                        'string. type "help pdf" for a list.'],form)
                end
                if isnan(test)
                    error('invalid parameters given for hyperprior form %s.',form)
                end
                hyperparamForms{end+1} = form;
                hyperparamParams{end+1} = params;
                skipNext = true;
            else
                error(['input #%i not recognized.  note that ' ...
                    'hyperparameterParams inputs must be cell matrices ' ...
                    'of scalar numeric values.'],nNamedInp+v+1)
            end
        elseif isnumeric(varargin{v}) && isscalar(varargin{v})
            if mod(varargin{v},1) == 0
                NS = varargin{v};
            else
                error('input #%i not recognized.',nNamedInp+v)
            end
        elseif iscell(varargin{v}) && isvector(varargin{v})
            hasFunc = cellfun(@(x) isa(x,'function_handle'),varargin{v});
            hasNone = cellfun(@isempty,varargin{v});
            if all(hasFunc | hasNone)
                functions = varargin{v};
            else
                error('input #%i not recognized.',nNamedInp+v)
            end
        else
            error('input #%i not recognized.',nNamedInp+v)
        end
    end
end



nHyperparams = length(hyperparamForms);
%if not enough HP ranges given, missing ones are a repeat of last given.
hpRange = repmat(xRange(end,:),[nHyperparams 1]);
for n = 2:size(xRange,1)
    hpRange(n-1,:) = xRange(n,:);
end
xRange = xRange(1,:);

%do any distributions have more than 4 parameters???
% hyperparamColors = makecolormap('darkblue','lightblue',nHyperparams);
% jointColor = mean(hyperparamColors,1);
hyperparamColors = getcolors('plum','indigo','blue','green');
% hyperparamColors = getcolors('lightplum','lightindigo','lightblue','lightgreen');
hyperparamColors = hyperparamColors(1:nHyperparams,:);
jointColor = mean(hyperparamColors,1);
jointColor = mean([jointColor; jointColor; 0.8*[1 1 1]; 1 1 1],1);
% jointColor = getcolors('darkblue');
priorColor = getcolors('black');
simPriorColor = 0.55*[1 1 1];

tfsm = 1.05;   %TitleFontSizeMultiplier
lfsm = 1;   %LabelFontSizeMultiplier
NX = 1000;  %points over domain

%standardize names
priorForm = internal.stats.getDistributionName(priorForm);
for hp = 1:nHyperparams
	hyperparamForms{hp} = ...
        internal.stats.getDistributionName(hyperparamForms{hp});
end

%functions
if exist('functions','var')
    if ~isequal(length(functions),nHyperparams)
        error(['if functions is given, then the number of functions' ...
            'must be equal to the number of hyperparameters.'])
    end
else
    functions = cell([nHyperparams 1]);
end

%% set up for prior testing plot
if noNewFig,    clf; f = gcf;
else,           f = figure('color',[1 1 1]);
end

if ~strcmp(f.WindowStyle,'docked')
    fpos = f.Position;
    f.Position = [fpos([1 2]) [600 625]*figScaling];
end

nRows = 5; nCols = 2;
titleStr = tiledlayout(nRows,nCols,'tilespacing','normal','padding','compact');

%% hyperprior distributions
nexttile(1,[2 1]); hold on
if nHyperparams==1, tileTitle = title('hyperprior');
else,               tileTitle = title('hyperpriors'); end
hpLabels = cell([1 nHyperparams]);
for hp = 1:nHyperparams
    xx = linspace(hpRange(hp,1),hpRange(hp,2),NX);
    hpHandle = plot(xx, ...
        pdf(hyperparamForms{hp},xx,hyperparamParams{hp}{:}), ...
        'linewidth',linePt*1.5,'color',hyperparamColors(hp,:));
    hpHandle.Color = [hpHandle.Color 0.75];
    g = repmat('%g,',[1 length(hyperparamParams{hp})]);
    hpLabels{hp} = sprintf(['[%i]  %s(' g(1:end-1) ')'], ...
        hp,hyperparamForms{hp},hyperparamParams{hp}{:});
end
%format
legend(hpLabels,'location','best','box','off')
set(gca,'fontsize',fontSz)
set(gca,'TitleFontSizeMultiplier',tfsm,'LabelFontSizeMultiplier',lfsm)
% tileTitle.HorizontalAlignment = 'left';
% tileTitle.Position = [0 tileTitle.Position(2:end)];

%% hyperparameters (sampled from hyperpriors)
%sample hyperparameter values from the hyperpriors as specified ...
hyperparamValues = NaN([NS nHyperparams]);
for hp = 1:nHyperparams
    if isempty(functions{hp})
        hyperparamValues(:,hp) = ...
            random(hyperparamForms{hp},hyperparamParams{hp}{:},[NS 1]);
    else
        f = functions{hp};
        hyperparamValues(:,hp) = ...
            f(random(hyperparamForms{hp},hyperparamParams{hp}{:},[NS 1]));
    end
end
%... and plot their (univariate or bivariate) distribution
switch nHyperparams
  case 1
    %univariate
    nexttile(2,[2 1]); hold on
    tileTitle = title('hyperparameter draws');
    histogram(hyperparamValues, ...
        'facecolor',hyperparamColors(1,:), ...
        'facealpha',0.35,'edgecolor','none')
    xlabel('[1]','fontweight','normal')
  case 2
    %bivariate
    nexttile(2,[2 1]); hold on
    tileTitle = title('hyperparameter draws');
    %temporarily mess with plot
    axunits = get(gca,'units'); %original settings
    set(gca,'units','points');  %temporarily change settings
    axpos = get(gca,'pos');
    set(gca,'units',axunits);   %revert to original settings
    %marker size
    mPtSize = (axpos(3)*0.025)^2; %marker width will be exactly 2.5% of 
                                  %the x-axis width.
    scatter(hyperparamValues(:,1),hyperparamValues(:,2), ...
        mPtSize,'markerfacecolor',jointColor, ...
        'markerfacealpha',0.35,'markeredgecolor','none')
    %format
    if isempty(functions{1}),   xLabelStr = '[1]';
    else,                       xLabelStr = 'function of [1]';
    end
    xlabel(xLabelStr,'fontweight','normal')
    if isempty(functions{1}),   yLabelStr = '[2]';
    else,                       yLabelStr = 'function of [2]';
    end
    ylabel(yLabelStr,'fontweight','normal')
    set(gca,'fontsize',fontSz,'box','on')
  otherwise
    %too many to plot.
end
if ismember(nHyperparams,[1 2])
    %format
    set(gca,'fontsize',fontSz,'box','on')
    set(gca,'TitleFontSizeMultiplier',tfsm,'LabelFontSizeMultiplier',lfsm)
%     tileTitle.HorizontalAlignment = 'left';
%     tileTitle.Position = [0 tileTitle.Position(2:end)];
end

%% simulated prior distributions (simulated using sampled hyperparameters)
nexttile(2*nCols+1,[nRows-2 nCols]); hold on
tileTitle = title('simulated prior distributions');
xx = linspace(xRange(1),xRange(2),NX);
h = gobjects([2 1]);

%N simulated priors
simPriorProb = zeros([NS length(xx)]);
for n = 1:NS
    %grab the nth set of hyperparameters and simulate a prior
    simHyperparam = num2cell(hyperparamValues(n,:));
    simPriorProb(n,:) = pdf(priorForm,xx,simHyperparam{:});
    %... and plot this simulation
    h(1) = plot(xx,simPriorProb(n,:), ...
        'color',simPriorColor,'linewidth',linePt/2);
    h(1).Color(4) = 0.25; %transparent line
end
%plot the average of the simulated priors
h(2) = plot(xx,mean(simPriorProb),'-', ...
    'color',priorColor,'linewidth',linePt*2);

%~*~ reference prior ~*~
if includeRef
    %yes
    refPrior = pdf(priorForm,xx,referenceParams{:});
    h(3) = plot(xx,refPrior, ...
        'color','k','linestyle','--','linewidth',linePt*2);

    gg = repmat('%g,',[1 length(referenceParams)]);
    refLabel = sprintf(['reference prior: %s(' gg(1:end-1) ')'], ...
                priorForm,referenceParams{:});
    %     sprintf('%i possible priors\n      from hyperpriors',NS), ...
    %     sprintf('average prior simulated\n      from hyperpriors'), ...
    hLabels = {sprintf('%i simulated %s priors',NS,priorForm), ...
               'average of simulated priors', ...
               refLabel};
else
    %none given
    hLabels = {sprintf('%i simulated %s priors',NS,priorForm), ...
               'average of simulated priors'};
end
legend(h,hLabels,'location','southoutside','box','off')

%format
set(gca,'fontsize',fontSz)
set(gca,'TitleFontSizeMultiplier',tfsm,'LabelFontSizeMultiplier',lfsm)
% tileTitle.HorizontalAlignment = 'left';
% tileTitle.Position = [0 tileTitle.Position(2:end)];

%% y-axis max
% %replace Inf with NaN
% simprob(simprob==Inf) = NaN; simprob(isnan(simprob)) = max(simprob(:));
% if includeRef
%     refPrior(refPrior==Inf) = NaN; refPrior(isnan(refPrior)) = max(refPrior(:));
% end
% %start with just above the average prior's max
% avgPriorMax = max(mean(simprob))*1.5
% ymax = avgPriorMax;
% if includeRef
%     refPriorMax = max(refPrior)*1.5
%     ymax = max(refPriorMax,avgPriorMax)
% end
% %get a good y-limit based on the simulated priors ...
% eachPriorMax = max(mean(simprob,2));
% for sc = 1:6
%     if mean(any(simprob>eachPriorMax*sc,2)) > 0.8
%         ymax = eachPriorMax*sc;
%         sc
%     end
% end
% %... and the reference prior
% if includeRef
%     if ymax > max(refPrior)
%         if ymax > max(refPrior)/6
%             ymax = ymax/2;
%         end
%     else
% %         for sc = 2:6
% %             if ymax < max(refPrior)*sc
% %                 ymax = max(refPrior)*sc;
% %                 'ref'
% %                 sc
% %             end
% %         end
%     end
% end
% ylim([0 ymax])

end