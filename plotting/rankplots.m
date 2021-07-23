function rankplots(samples,varargin)
%RANKPLOTS displays and a marginal density for each parameter.
% 
% this function will generate a single figure with one subplot per chain,
% where is subplot is a histogram of that chain's posterior samples' rank
% within the pool of samples across all chains.
% 
% if all chains are behaving similarly, the rank plot for each and every
% chain should appear uniform. 
% 
% this new diagnostic plot was suggested in:
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% 
% a maximum of 25 figures may be generated in a single call, to minimize
% the risk of memory overload.   
% 
% 
% RANKPLOTS(SAMPLES)
% RANKPLOTS(SAMPLES,PARAMETERREQUEST)
%   SAMPLES is the structure of posterior samples (in the format generated
%   by extractsamples.m).  
%   
%   if SAMPLES is the only input, then a figure for all parameters/
%   parameter instances will (attempt to) be generated.  
%   if PARAMETERREQUEST is also given, then a figure for only the requested
%   parameters/parameter instances will be generated.   
%   
%   if PARAMETERREQUEST is a parameter name string, then a figure will
%   be generated for [all instances of] that parameter only. 
%   if PARAMETERREQUEST is a cell of parameter name strings, then a
%   traceplot will be generated  for [all instances of] each of those
%   parameters. 
%   
%   PARAMETERREQUEST may also include particular parameter instances, by 
%   including a valid index or comma-separated indices in brackets 
%   after the parameter name (e.g., 'mu[3]','sigma[2,1]').
%   
%   * may also be used as a wildcard to include any & all matching
%   parameter instance names (e.g., '*_alpha', 'beta[*,1]').
% 
% 
% RANKPLOTS(...,NBINS)
%   a custom number of bins, NBINS, may be specified.
% 
% 
% See also TRACEDENSITY
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

matstanlib_options

%% parse inputs
%samples
if ~isstruct(samples)
    error(['first input must be a structure of posterior samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end

%nBins
if isempty(varargin) || ~any(cellfun(@isnumeric,varargin))
    nBins = []; %determine # of bins from # of samples
else
    isNumeric = find(cellfun(@isnumeric,varargin));
    if length(isNumeric)==1
        nBins = varargin{isNumeric};
        if isscalar(nBins) && ~(mod(nBins,1)>0) && nBins>=10
            varargin(isNumeric) = []; %pop out of varargin
        else
            error('nBins input must be an integer >= 10.')
        end
    elseif length(isNumeric)>1
        error('can only have one numeric-type input (nBins).')
    else
        %no nBins given
    end
end

%parameterRequest
if isempty(varargin) || ...
        ~(isempty(varargin{1}) || ischar(varargin{1}) || iscell(varargin{1}))
    parameterRequest = fieldnames(samples); %default
else                                  %pop out of varargin
    parameterRequest = varargin{1};   varargin = varargin(2:end);
    if isempty(parameterRequest)
        parameterRequest = fieldnames(samples); %default
    elseif ischar(parameterRequest)
        parameterRequest = {parameterRequest};
    elseif iscell(parameterRequest) 
        if ~all(cellfun(@ischar,parameterRequest))
            error(['parameterRequest must be input as a string or ' ...
                   'cell of strings.'])
        end
    end
end

%
if ~isempty(varargin)
    error('too many inputs.')
end

%% prepare to generate plots
chainColors = getcolors('mcmc');
barColor = getcolors('lightgray');

%maximum number of figures to generate
maxNfigures = 25;

%create a list of parameter instances
parameters = getparaminstances(parameterRequest, ...
    fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));

%warn if over the figure limit
if length(parameters) > maxNfigures
    warning(['requested %i traceplots, but this function will only ' ...
        'generate up to %i traceplots in a single call (to protect ' ...
        'against memory overload and subsequent freezing). ' ...
        'try requesting fewer parameters/parameter instances.  ' ...
        'continuing to generate the first %i traceplots now ...'], ...
        length(parameters),maxNfigures,maxNfigures)
end

%% make one figure for each parameter, with one rank plot per chain
isInstance = cellfun(@(x) any(x=='['),parameters);
for p = 1:min([length(parameters),maxNfigures])
    %account for parameters vs. parameter instances
    if isInstance(p)
        [parameter,ind] = str2ind(parameters{p});
        chains = samples.(parameter)(:,:,ind{:});
    else
        parameter = parameters{p};
        chains = samples.(parameter);
    end
    
    %extract info
    [nIterations,nChains] = size(chains,[1 2]);
    nSamples = nIterations*nChains;
    
    %rank (relative to **all** samples)
    sz = size(chains);  %extract size
    chains = tiedrank(chains(:));    %rank
    chains = reshape(chains,sz);  %convert from back to matrix
    
    %start a figure ...
    dumf = figure(999); %dummy figure to protect sizing
    f = figure('color',[1 1 1]);
    fpos = f.Position;
    close(dumf.Number); %close dummy figure
    %... and a layout
    t = tiledlayout('flow','tilespacing','compact');
    t.Title.String = parameters{p};
    set(t.Title,'fontweight','bold','interpreter','none','fontsize',fontSz*1.1)
    
    %%% rank plot (one per chain) %%%
    for m = 1:nChains
        ax = nexttile; hold on;
        if m==1
            if ~isempty(nBins)
                binOpts = {nBins}; %(use input)
            elseif nSamples < 500
                binOpts = {10}; %minimum of 10 bins
            elseif nSamples < 1250
                binOpts = {20}; %use 20 bins if between 500 and 1249 samples
            elseif nSamples < 2500
                binOpts = {35}; %use 35 bins if between 1250 and 2499 samples
            else
                binOpts = {50}; %use 50 bins if at least 2500 samples
            end
        end
        histogram(chains(:,m),binOpts{:}, ...
            'facecolor',chainColors(m,:),'edgecolor','k')
        %format
        xlim([-0.1*nSamples 1.1*nSamples]); xt = get(gca,'xtick'); 
        xt(xt<0) = []; xt(xt>nSamples) = []; set(gca,'xtick',xt)
        ax.YAxis.Visible = 'off'; %y-axis is uninformative
        % xlabel('rank')
        % xlabel(sprintf('chain %i',m))
        title(sprintf('chain %i',m),'fontweight','normal')
        set(gca,'fontsize',fontSz)
    end
    
    %resize the figure based on the number of tiles
    gridSize = t.GridSize;
    f.Position = [fpos(1:2) (gridSize.*[300 150] + [0 75])*figScaling];
end    

end
