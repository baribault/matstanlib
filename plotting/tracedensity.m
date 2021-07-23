function tracedensity(samples,varargin)
%TRACEDENSITY displays chain traces and a marginal density for each parameter.
% 
% this function will generate one figure per parameter, where each figure
% contains a chain trace plot and a marginal density histogram.  
% optionally, diagnostic overlays may be included. 
% 
% a maximum of 25 trace plots may be generated in a single call, to
% minimize the risk of memory overload.  
% 
% 
% TRACEDENSITY(SAMPLES)
% TRACEDENSITY(SAMPLES,PARAMETERREQUEST)
%   SAMPLES is the structure of posterior samples (in the format generated
%   by extractsamples.m).  
%   
%   if SAMPLES is the only input, then a figure for all parameters/
%   parameter instances will (attempt to) be generated.  
%   if PARAMETERREQUEST is also given, then a figure for only the requested
%   parameters/parameter instances will be generated.   
%   
%   if PARAMETERREQUEST is a parameter name string, then a traceplot will
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
% TRACEDENSITY(...,'ranknorm',...)
%   adding the string 'ranknorm' to the inputs 
% 
% 
% TRACEDENSITY(..., DIAGNOSTICS)
%   if the DIAGNOSTICS structure is also given, rug plots of divergent
%   transitions (if any occured) and a rugplot of iterations where the
%   maximum tree depth was hit (if any), will be added to the bottom of the
%   plot.  
% 
% TRACEDENSITY(..., DIAGNOSTICS,MAXTREEDEPTH)
%   if a maximum tree depth other than the default value of 10 was used
%   during sampling with the NUTS algorithm, it should be given as a final
%   optional input, MAXTREEDEPTH. 
% 
% 
% See also RANKPLOTS, EXTRACTSAMPLES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

matstanlib_options

%% parse inputs
%samples
if ~isstruct(samples)
    error(['first input must be a structure of posterior samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end

%rankNormalize
rankNormalize = false;
if ~isempty(varargin)
    isRankNormStr = ...
        cellfun(@(x) ischar(x) && strcmp(x,'ranknorm'),varargin);
    if any(isRankNormStr)
        rankNormalize = true;
        varargin(isRankNormStr) = []; %pop out of varargin
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

overlayDx = false;
if ~isempty(varargin)
    %diagnostics
    if isstruct(varargin{1})
        if all(isfield(varargin{1},{'divergent__','treedepth__'}))
            overlayDx = true;
            divergent = any(varargin{1}.divergent__,2);
            treedepth = max(varargin{1}.treedepth__,[],2);
            maxTD = 10; %default
        else
            error(['an optional struct-type input was given, but it ' ...
            'is missing expected fields for the diagnostics structure.'])
        end
    else
        error(['an optional input was given, but it was not struct-type, ' ...
            'when the diagnostics structure was expected.'])
    end
    if length(varargin)==1 && any(treedepth(:) > 10)
        error(['treedepth value above the default value of 10 were found, ' ...
            'but the custom maximum tree depth setting was not given as ' ...
            'the final optional input.'])
    end
    if length(varargin) > 1
        if ~isnumeric(varargin{2}) || ~isscalar(varargin{2}) || ...
                mod(varargin{2},1) > 0 && varargin{2} > 1
            error('max tree depth must be a single integer > 1.')
        else
            maxTD = varargin{2};
        end
    end
    if length(varargin) > 1
        error('too many inputs.')
    end
end
clearvars varargin

%% prepare to generate plots
chainColors = getcolors('mcmc');
divergentColor = getcolors('red');
treedepthColor = getcolors('orange');

%maximum number of figures to generate
maxNfigures = 25;

%determine positioning within each figure
% ax1pos = [0.075 0.235 0.7  0.635]; %trace plot
% ax2pos = [0.8   0.235 0.12 0.635]; %marginal density
ax1pos = [0.1   0.235 0.75 0.7]; %trace plot
ax2pos = [0.875 0.235 0.1  0.7]; %marginal density

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

%prepare for diagnostic overlays
if overlayDx
    %convert treedepth to an indicator of hitting max tree depth
    treedepth = treedepth==maxTD;
    %which diagnostics to plot?
    quantities = [1 2];
    if ~any(divergent(:)), quantities(1) = 0; end
    if ~any(treedepth(:)), quantities(2) = 0; end
    quantities(quantities==0) = [];
    if isempty(quantities)
        %nothing to plot after all!
        overlayDx = false;
    end
end

%% make trace plots for each parameter
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
    
    %rank normalize?
    if rankNormalize
        chains = ranknorm(chains);
    end
    
    %extract info
    nIterations = size(chains,1);
    nChains = size(chains,2);
    postMean = mean(chains(:));
    
    %start a figure ...
    dumf = figure(999); %dummy figure to protect sizing
    f = figure('color',[1 1 1]);
    fpos = f.Position;
    f.Position = [fpos(1:2) [600 380]*figScaling];
    close(dumf.Number); %close dummy figure
    
    %%% traceplot (plot chains) %%%
    ax1 = axes('position',ax1pos);
    hold on
    %dummy thick lines for legend
    h = gobjects([1 nChains]);
    for n = 1:nChains
        h(n) = plot(-1,chains(1,1),'color',chainColors(n,:),'linewidth',linePt*3);
    end
    %actual traces
    for n = 1:size(chains,2)
        plot(chains(:,n),'color',chainColors(n,:))%,'linewidth',linePt)
    end
%     %title, with Rhat
%     rhat = splitrhat(chains);
%     protectedParam = strrep(parameters{p},'_','\_'); 
%     title({['\bf ' protectedParam],['\rm ' sprintf('Rhat = %.3f',rhat)]}, ...
%         'fontweight','normal','interpreter','tex'), with Rhat
    protectedParam = strrep(parameters{p},'_','\_'); 
    title(protectedParam,'interpreter','tex')
    %format x-axis
    xlim([0 nIterations+1])
    xlabel('iteration')
    
    %%% rug plots (plot divergences, hit max treedepth) %%%
    if overlayDx
        %expand the y-axis to fit the diagnostics
        xl = xlim;
        yl = ylim;
        barMarkerSz = markSz*0.75;
        rug_pos(1) = min(chains(:)) - 0.01*barMarkerSz*diff(yl);
        rug_pos(2) = min(chains(:)) - 0.02*barMarkerSz*diff(yl);
        rug_pos(3) = min(chains(:)) - 0.03*barMarkerSz*diff(yl);
        %make one rugplot per diagnostic
        c = 0;
        for q = quantities
            if     q==1, dx = divergent; dxColor = divergentColor;
            elseif q==2, dx = treedepth; dxColor = treedepthColor;
            end
            %rug plot
            c = c + 1; %increment counter
            x_dx = find(dx);
            y_dx = rug_pos(c)*ones(size(x_dx));
            plot(x_dx,y_dx,'|','markersize',barMarkerSz,'color',dxColor)
        end
        %reestablish limits, just to be sure
        xlim(xl); ylim([mean(rug_pos(c:c+1)) yl(2)]);
    end
    %format
    legend(h,arrayfun(@num2str,1:nChains,'uni',0), ...
        'location','southoutside','orientation','horizontal')
    set(ax1,'pos',ax1pos) %reposition!
    set(ax1,'fontsize',fontSz,'box','on')
    
    
    %%% marginal distribution %%%
    ax2 = axes('position',ax2pos);
    hold on
    %histogram of samples
    histogram(chains(:),'orientation','horizontal', ...
        'facecolor',getcolors('gray'),'edgecolor',getcolors('gray'))
    %line at posterior mean
    plot(xlim,postMean*[1 1],'k-','linewidth',linePt)
    %line at posterior median
    plot(xlim,postMean*[1 1],'k--','linewidth',linePt)
    %format
%     set(ax2,'ylim',ax1.YLim,'ytick',ax1.YTick);
    set(ax2,'ylim',ax1.YLim,'ytick',ax1.YTick,'yticklabels',{});
    set(ax2,'xtick',[],'box','off')
    set(ax2,'fontsize',fontSz,'xcolor','w')
%     title({'mean =',sprintf('%.3g',postMean)}, ...
%         'fontweight','normal','fontsize',fontSz-2)
    xlabel({'marginal','posterior','density'},'color','k', ...
        'fontsize',fontSz-2)
end    

end
