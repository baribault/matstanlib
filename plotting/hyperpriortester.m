function hyperpriortester(xRange,priorForm,idealPriorParams,varargin)
% HYPERPRIORTESTER helps visualize hyperpriors' effects on a prior.
% 
% *** REMEMBER: MATLAB AND STAN HAVE DIFFERENT PARAMETERIZATIONS ***
% ****          FOR SOME DISTRIBUTIONS!!!          (e.g., Gamma) ***
% 
% this function will test the effect of different hyperprior distributions
% on possible parameter values at the prior level, by visualizing how
% different hyperpriors distributions combine to constrain the prior. 
% 
% specifically, this function will generate a figure your candidate
% hyperprior distribution (with the given parameterization), simulate thier
% joint distribution (to foresee challenges), and will contrast your ideal
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
% HYPERPRIORTESTER(XRANGE,PRIORFORM,IDEALPRIORPARAMS, ...
%                  HYPERPRIORFORM,HYPERPRIORPARMS, ...
%                  [HYPERPRIORFORM2,HYPERPRIORPARMS2] [, ...])
%   XRANGE is the range of parameter values to consider (and plot!). 
%   the same values will be used for the hyperparameter plots, unless
%   otherwise given with additional rows in XRANGE. 
%   
%   PRIORFORM indicates what form of the prior you are using is (i.e., the
%   distribution name).   
%   
%   IDEALPRIORPARAMS contains the parameters that would define the ideal
%   prior distribution you'd use (i.e., what hyperparameters you'd use if
%   you were to use this prior non-hierarchically).  
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
%   each hyperparameter before using it for simulation.  
%   (if there is no function to apply
%   these functions can add a minimum value, scale, etc.  the only
%   requirement is that the functions only have ONE input (the parameter
%   value!).
%   NOTE: these functions are only used to get the expected resultant
%   prior!
% 
% 
% Example:
%   > f = {@(x) x + 1; @(x) x}; %adds 1 to shape parameter
%   > hyperpriortester([0 10],'gamma',{5,1/1}, ... %ideal prior
%                             'gamma',{5,1/2}, ... %hyperprior for shape
%                             'gamma',{2,1/1},f)   %hyperprior for rate/scale
% 
% 
% See also PDF
% 
% (c) beth baribault 2020 ---                                 > matstanlib 

matstanlib_options

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
if ~iscell(idealPriorParams) || ~isvector(idealPriorParams)
    error('idealPriorParams must be a cell type vector.')
elseif ~all(cellfun(@isnumeric,idealPriorParams) & cellfun(@isscalar,idealPriorParams))
    error('idealPriorParams must be a cell type vector of scalar numeric values.')
end

%prior form AGAIN
if ~ischar(priorForm)
    error('priorForm must be a string.')
else
    try
        test = random(priorForm,idealPriorParams{:});
    catch
        error(['priorForm must be a valid distribution string.  ' ...
            'type "help pdf" for a list.'])
    end
    if isnan(test)
        error('invalid parameters given for prior form %s.',priorForm)
    end
end

%% parse optional inputs

%initialize lists
hyperparamForms = {};
hyperparamParams = {};

%defaults
noNewFig = false;
parameterName = '';
NX = 1000; %cannot be changed
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
        elseif ischar(varargin{v})
            parameterName = varargin{v};
        elseif isnumeric(varargin{v}) && isscalar(varargin{v})
            if mod(varargin{v},1) == 0
                NX = varargin{v};
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
hyperparamColors = getcolors('blue','green','yellow','mandarin');
hyperparamColors = hyperparamColors(1:nHyperparams,:);
jointColor = getcolors('darkblue');
priorColor = getcolors('black');
simPriorColor = 0.55*[1 1 1];

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
    if nHyperparams == 1
        f.Position = [fpos([1 2]) [900 250]*figScaling];
    else
        f.Position = [fpos([1 2]) [750 520]*figScaling];
    end
end
%subplots, ugh
if nHyperparams == 1
    sp1 = 1; sp2 = 3; plt3 = [2 3]; loc = 'eastoutside';
else
    sp1 = 2; sp2 = 2; plt3 = [3 4]; loc = 'eastoutside';
end

%%% hyperpriors, curves
subplot(sp1,sp2,1); hold on
title('hyperpriors')
hpLabels = cell([1 nHyperparams]);
for hp = 1:nHyperparams
    xx = linspace(hpRange(hp,1),hpRange(hp,2),NX);
    plot(xx,pdf(hyperparamForms{hp},xx,hyperparamParams{hp}{:}), ...
        'linewidth',linePt,'color',hyperparamColors(hp,:))
    g = repmat('%g,',[1 length(hyperparamParams{hp})]);
    hpLabels{hp} = sprintf(['[%i]  %s(' g(1:end-1) ')'], ...
        hp,hyperparamForms{hp},hyperparamParams{hp}{:});
end
legend(hpLabels,'location','best','box','off')
set(gca,'fontsize',fontSz)

%%% hyperpriors, scatter
if nHyperparams > 1
    subplot(2,2,2); hold on
    %temporarily mess with plot
    axunits = get(gca,'units'); %original settings
    set(gca,'units','points');  %temporarily change settings
    axpos = get(gca,'pos');
    set(gca,'units',axunits);   %revert to original settings
    %marker size
    mPtSize = (axpos(3)*0.025)^2; %marker width will be exactly 2.5% of 
                                  %the x-axis width.
    scatter( ...
        random(hyperparamForms{1},hyperparamParams{1}{:},[NX 1]), ...
        random(hyperparamForms{2},hyperparamParams{2}{:},[NX 1]), ...
        mPtSize,'markerfacecolor',jointColor, ...
        'markerfacealpha',0.35,'markeredgecolor','none')
    xlabel('[1]','fontweight','bold')
    ylabel('[2]','fontweight','bold')
    set(gca,'fontsize',fontSz,'box','on')
end

%%% prior
subplot(sp1,sp2,plt3); hold on
t = 'prior'; if ~isempty(parameterName), t = [t '(' parameterName ')']; end
title(t,'interpreter','none')
xx = linspace(xRange(1),xRange(2),NX);
h = gobjects([2 1]);

%simulated expected prior ...
simprob = zeros([NS length(xx)]);
for n = 1:NS
    fromHyperprior = cell([nHyperparams 1]);
    for hp = 1:nHyperparams
        if isempty(functions{hp})
            fromHyperprior{hp} = ...
                random(hyperparamForms{hp},hyperparamParams{hp}{:});
        else
            f = functions{hp};
            fromHyperprior{hp} = ...
                f(random(hyperparamForms{hp},hyperparamParams{hp}{:}));
        end
    end
    simprob(n,:) = pdf(priorForm,xx,fromHyperprior{:});
    %... and plot all simulations
    h(1) = plot(xx,simprob(n,:), ...
        'color',simPriorColor,'linewidth',linePt/2);
    h(1).Color(4) = 0.25; %transparent line
end
% ... and the average simulation
h(2) = plot(xx,mean(simprob),'-', ...
    'color',priorColor,'linewidth',linePt*2);

%~*~ ideal prior ~*~
idealPrior = pdf(priorForm,xx,idealPriorParams{:});
h(3) = plot(xx,idealPrior, ...
    'color','k','linestyle','--','linewidth',linePt*2);

gg = repmat('%g,',[1 length(idealPriorParams)]);
idealLabel = sprintf(['ideal prior: %s(' gg(1:end-1) ')'], ...
            priorForm,idealPriorParams{:});
%     sprintf('%i possible priors\n      from hyperpriors',NS), ...
%     sprintf('average prior simulated\n      from hyperpriors'), ...
legend(h, ...
     sprintf('%i simulated priors',NS),'average of simulated priors',idealLabel, ...
    'location',loc,'box','off')
set(gca,'fontsize',fontSz)

smax = max(mean(simprob));
if smax > max(idealPrior)*2.5
    ymax = max(idealPrior)*2.5;
else
    ymax = max(idealPrior)*1.35;
end
ylim([0 ymax])

end