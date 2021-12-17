function rankplots(samples,varargin)
%RANKPLOTS displays and a marginal density for each parameter.
% 
% this function will generate a single figure with one subplot per chain,
% where each subplot is a histogram of that chain's posterior samples' rank
% within the pool of samples across all chains.  
% 
% if all chains are behaving similarly, the rank plot for each and every
% chain should appear uniform.  the faded black line represents a perfectly
% uniform distribution.  
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
% See also TRACEDENSITY, PLOTESS, PLOTMCSE
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

msl.options

%% parse inputs
%samples
if isstruct(samples)
    %gave samples
elseif isnumeric(samples) && ismatrix(samples)
    %gave chains
    chainParam = '';
    if ~isempty(varargin)
        isParamStr = ...
            cellfun(@(x) ischar(x) && ~strcmp(x,'ranknorm'),varargin);
        if any(isParamStr)
            chainParamInd = find(isParamStr,1);
            if ~isempty(varargin{chainParamInd})
                chainParam = varargin{chainParamInd};
            end
            varargin(isParamStr) = []; %pop out of varargin
        end
    end
    if isempty(chainParam), chainParam = 'NO_TITLE'; end
    %convert to struct
    tmp_samples.(chainParam) = samples;
    clearvars samples
    samples = tmp_samples;
else
    error(['first input must be a structure of posterior samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end
% % % %samples
% % % if ~isstruct(samples)
% % %     error(['first input must be a structure of posterior samples ' ...
% % %         '(consisitent with the output of extractsamples.m).'])
% % % end

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
    if p==1
        [nIterations,nChains] = size(chains,[1 2]);
        nSamples = nIterations*nChains;
        %if not given, determine number of bins to use in histogram
        if isempty(nBins)
            if nSamples <= 500
                nBins = 10; %minimum of 10 bins
            elseif nSamples <= 2000
                nBins = 20; %use 20 bins if between 501 and 2000 samples
            elseif nSamples <= 8000
                nBins = 40; %use 40 bins if between 2001 and 8000 samples
            else
                nBins = 50; %use 50 bins if more than 8000 samples
            end
        end
    end
    
    %rank (relative to **all** samples)
    sz = size(chains);              %extract size
    chains = tiedrank(chains(:));   %rank
    chains = reshape(chains,sz);    %convert from back to matrix
    
    %start a figure ...
    dumf = figure(999); %dummy figure to protect sizing
    f = figure('color',[1 1 1]);
    fpos = f.Position;
    f.Position = [fpos(1:2) [520 420]*figScaling];
    close(dumf.Number); %close dummy figure
    %... and a layout
    t = tiledlayout('flow','tilespacing','compact','padding','compact');
    t.Title.String = parameters{p};
    set(t.Title,'fontweight','bold','interpreter','none','fontsize',fontSz*1.1)
    
    %%% rank plot (one per chain) %%%
    ax = gobjects([nChains 1]);
    for m = 1:nChains
        ax(m) = nexttile; hold on;
        %underlay line to represent uniformly distributed ranks
        unifLineOpts = {'color','k','linewidth',linePt};
        hl0 = plot([-0.1*nSamples -0.025*nSamples],[1 1]*nIterations/nBins, ...
            unifLineOpts{:},'linestyle',':');
        hl = plot([-0.025*nSamples 1.025*nSamples],[1 1]*nIterations/nBins, ...
            unifLineOpts{:},'linestyle','-');
        hl1 = plot([1.025*nSamples 1.1*nSamples],[1 1]*nIterations/nBins, ...
            unifLineOpts{:},'linestyle',':');
        hl0.Color = [hl0.Color 0.5];
        hl.Color = [hl.Color 0.5];
        hl1.Color = [hl1.Color 0.5];
        %rank plot
        binOpts = {nBins};
        histogram(chains(:,m),binOpts{:}, ...
            'facecolor',chainColors(m,:),'edgecolor','k', ...
            'facealpha',1)
        %format
        xlim([-0.1*nSamples 1.1*nSamples]); 
        ax(m).YAxis.Visible = 'off'; %y-axis is uninformative
        % xlabel('rank')
        % xlabel(sprintf('chain %i',m))
        title(sprintf('chain %i',m),'fontweight','normal')
        set(gca,'fontsize',fontSz)
        %LAST ORDER OF BUSINESS is xticks; otherwise yucky angles
        if m==1
            xt = get(gca,'xtick');
            xt(xt<0 | xt>nSamples) = [];
        end    
        set(gca,'xtick',xt)          %x-axis ticks are same for all
    end
    linkaxes(ax,'x') %zooming one zooms all
    
    %resize the figure based on the number of tiles
    gridSize = t.GridSize;
    f.Position = [fpos(1:2) (gridSize.*[300 150] + [0 75])*figScaling];
end

end
