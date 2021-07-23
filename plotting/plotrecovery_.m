function varargout = plotrecovery(estimatedValues,trueValues, ...
    parameterRequest,varargin)
%PLOTRECOVERY plots true parameter values against estimated values.
%   
% this function may be used to visualize the quality of parameter 
% recovery when a model is applied to simualted data.  by default, a 
% single figure is generated, and the parameters to be included on 
% the plot are automatically determined.
% 
% the plot will include circular markers representing point estimates, 
% vertical error bars representing central 95% credible intervals, 
% and an overlaid diagonal line representing perfect recovery. 
% 
% NOTE: this function assumes that data was simulated from known 
% parameter values (i.e., true values), and that getsamplestats.m has 
% already been called to create the estimatedValues structure. 
%   
% PLOTRECOVERY(ESTIMATEDVALUES,TRUEVALUES)
%   generates a recovery plot for all parameters based on ESTIMATEDVALUES, 
%   a structure of estimated values, and TRUEVALUES, a structure of true
%   parameter values used to simulate data.  
% 
% 
% RECOVERYCOUNTS = PLOTRECOVERY(ESTIMATEDVALUES,TRUEVALUES,PARAMETERREQUEST)
%   PARAMETERREQUEST and RECOVERYCOUNTS are an optional input and output.
%   
%   if PARAMETERREQUEST is given, then the recovery plot will be generated
%   for the requested parameters *only*. PARAMETERREQUEST must be a valid
%   parameter name string or a cell of valid parameter name strings.
%   if PARAMETERREQUEST is an empty cell, all parameters will be plotted
%   (the default behavior). 
%   
%   optionally, this function can return RECOVERYCOUNTS, a two-element
%   vector including the number of parameters successfully and
%   unsuccessfully recovered, respectively. 
%   note that a parameter instance is considered "successfully recovered"
%   if the true parameter value is within the estimated 95% credible
%   interval.   
%   as such, RECOVERYCOUNTS may be used to determine the overall % of
%   parameters successfully recovered --- which should be approximately 95%
%   across all model parameters. 
% 
% 
% PLOTRECOVERY(...,POINTESTTYPE,MAKEERRORBARS,CREDINTTYPE,HIGHLIGHTVALUE,MARKERCOLOR)
%   four additional optional inputs may be given in any order (as long as
%   they are given after PARAMETERREQUEST). 
% 
%   POINTESTTYPE is a string indicating what posterior statitistic should
%   be used as the point parameter estimate.  valid options include:
%           'mean'    >>  mean of the posterior samples
%           'median'  >>  median of the posterior samples
%           'mode'    >>  mode of the posterior samples
%           'MAP'     >>  maximum a posteriori
%   note that these most already be fields in ESTIMATEDVALUES.
%   (by default, the 'mean' field in ESTIMATEDVALUES is used.) 
%   
%   MAKEERRORBARS must be a single logical value. if MAKEERRORBARS is true,
%   95% credible interval error bars are included. (by default, error bars
%   are omitted.)  
%   
%   CREDINTTYPE is a string indicating what type of CIs to plot (if any).
%   if MAKEERRORBARS is true, then the type of credible interval bars is
%   determined by CREDINTTYPE:
%           'CI'    >>  central 95% credible intervals
%           'HDI'   >>  95% highest density intervals
%   (by default, the credible interval type is central 95% intervals.)  
%   
%   HIGHLIGHTVALUE is a single numeric value.  
%   if HIGHLIGHTVALUE is given, a vertical line & a horizontal line will be
%   overlaid at the requested value. if HIGHLIGHTVALUE is empty, no value
%   is highlighted. (by default, no value is highlighted.) 
%   
%   MARKERCOLOR can offer a [nParameters 3] matrix of RGB-01 colors, one
%   for each recovered parameter. 
% 
% See also GETSAMPLESTATS, COLLECTTRUEVALUES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

matstanlib_options

%% example usage

% recoveryCounts = [0 0]
% recoveryCounts = recoveryCounts + ...
%           plotrecovery(estimatedValues,trueValues,{'parameter1of3'});
% recoveryCounts = recoveryCounts + ...
%           plotrecovery(estimatedValues,trueValues,{'parameter2of3'});
% recoveryCounts = recoveryCounts + ...
%           plotrecovery(estimatedValues,trueValues,{'parameter3of3'});

% OR

% recoveryCounts = plotrecovery(estimatedValues,trueValues); %all parameters

% THEN

% %recovery report
% proportionRecovered = recoveryCounts(1)/sum(recoveryCounts);
% fprintf('%i of %i parameters (%.2f%%) were sucessfully recovered!\n', ...
%     recoveryCounts(1),sum(recoveryCounts),proportionRecovered*100)

%% parse inputs
if ~isstruct(estimatedValues) || ~isstruct(trueValues)
    error(['first and second input must be the estimatedValues and ' ...
        'trueValues structures, respectively.'])
end

%parameterRequest: request a subset of parameters to plot
allValidParameters = intersect(fieldnames(estimatedValues),fieldnames(trueValues));
if nargin < 3 || isempty(parameterRequest)
    parameterRequest = allValidParameters;
elseif ischar(parameterRequest)
    parameterRequest = {parameterRequest};
elseif iscell(parameterRequest) && all(cellfun(@ischar,parameterRequest))
    %all good
else
    error(['the third input must be parameterRequest.  the requested ' ...
        'parameter name(s) must be input as a string or cell of strings ' ...
        '(or an empty cell, to request all parameters).'])
end
if any(~ismember(parameterRequest,allValidParameters))
    missingRequests = setdiff(parameterRequest,allValidParameters);
    error(['the following elements of parameterRequest were ' ...
        'not fields in the estimatedValues structure and/or ' ...
        'were not fields in the trueValues structure: %s'], ...
        ['''' strjoin(missingRequests,''', ''') ''''])
end

%% parse optional inputs
validPointEstTypes = {'mean','median','mode','MAP'};
pointEstTypeFound = false;      pointEstType = 'mean'; %default
validCredIntTypes = {'CI','CI95','CI50','HDI'};
credIntTypeFound = false;       credIntType = 'CI';    %default
makeErrorBarsFound = false;     makeErrorBars = false; %default
highlightValueFound = false;    highlightValue = [];   %default
minBoundsFound = false;
markerColorFound = false;

%loop over optional inputs
for v = 1:length(varargin)
    %pointEstType
    if ischar(varargin{v})
        if ismember(varargin{v},validPointEstTypes)
            if pointEstTypeFound
                error('can only have ONE optional pointEstType input string.')
            else
                pointEstType = varargin{v};
                pointEstTypeFound = true;
            end
        elseif ismember(upper(varargin{v}),validCredIntTypes) %case insensitive
            if credIntTypeFound
                error('can only have ONE optional pointEstType input string.')
            else
                credIntType = upper(varargin{v}); %case insensitive
                credIntTypeFound = true;
            end
        else
            error('unrecognized optional string type input.')
        end
    %makeErrorBars
    elseif islogical(varargin{v})
        if makeErrorBarsFound
            error(['can only have ONE optional logical input (makeErrorBars: ' ...
                'a logical value indicating whether or not to plot ' ...
                'credible intervals.)'])
        else
            makeErrorBars = varargin{v};
            makeErrorBarsFound = true;
        end
    % ... 
    elseif isnumeric(varargin{v})
        if isempty(varargin{v})
            %skip this input
    %highlightValue
        elseif isscalar(varargin{v})
            if highlightValueFound
                error('can only have ONE optional scalar numeric input (highlightValue).')
            else
                highlightValue = varargin{v};
                highlightValueFound = true;
            end
    %minBounds
        elseif numel(varargin{v})==2
            if diff(varargin{v}) <= 0
                error('minBounds must contain two *ascending* numeric elements.')
            elseif minBoundsFound
                error('can only have ONE optional two-element numeric input (minBounds).')
            else
                minBounds = varargin{v};
                minBoundsFound = true;
            end
    %markerColor: a vector or matrix of colors
        elseif ((isvector(varargin{v}) && length(varargin{v})==3) || ...
                (ismatrix(varargin{v}) && size(varargin{v},2)==3)) && ...
                all(varargin{v}(:) >= 0 & varargin{v}(:) <= 1)
            if markerColorFound
                error('can only have ONE optional input offering a [matrix of] colors.')
            else 
                if iscolumn(varargin{v}), varargin{v} = varargin{v}'; end
                if isrow(varargin{v}) || size(varargin{v},1)==length(parameterRequest)
                    markerColor = varargin{v};
                    markerColorFound = true;
                else
                    error(['markerColor input was %ix3 but should be %ix3 ' ...
                        '(for %i parameters).'],size(varargin{v},1), ...
                        length(parameterRequest),length(parameterRequest))
                end
            end
    %unrecognized (numeric) type!
        else
            error(['optional numeric input was given, but was not a ' ...
                'scalar (highlightValue) or matrix of colors (markerColor).'])
        end
    %unrecognized type!
    else
        error(['at least one optional input was given ' ...
            'but was not a valid, recognized type.'])
    end
end

%% prepare to plot
% ... just rename at this point
parameters = parameterRequest;

%get some colors
recoveredColor = getcolors('gray');    %cred int bars
failedColor = getcolors('red');        %cred int bars
if exist('markerColor','var')
    if isvector(markerColor)
        %same color for all parameters
        paramColors = repmat(markerColor,[length(parameters) 1]);
    else
        %different color for each parameter
        paramColors = markerColor;
    end
else
    paramColors = getcolors('ocean');           %point estimate markers
    if length(parameters) >= size(paramColors,1)
        paramColors = getcolors('mcmc');
        if length(parameters) >= size(paramColors,1)
            %so many parameters omg why
            paramColors = rand([length(parameters) 3]);
        end
    end
end

%and pick your CI type
switch credIntType
    case {'CI','CI95'}
        %95% central interval
        lowerBound = 'lowerCI';
        upperBound = 'upperCI';
    case 'CI50'
        %50% central interval
        lowerBound = 'lowerCI50';
        upperBound = 'upperCI50';
    case 'HDI'
        %95% highest denisty interval
        lowerBound = 'lowerHDI';
        upperBound = 'upperHDI';
end

%% recovery plot
%start a figure
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) [425 425]*figScaling];
close(dumf.Number); %close dummy figure
%start an axis
ax = axes;
hold on

%marker size ... sigh
%temporarily mess with plot
axunits = get(gca,'units'); %original settings
set(gca,'units','points');  %temporarily change settings
axpos = get(gca,'pos');
set(gca,'units',axunits);   %revert to original settings

%plot settings
mPtSize = (axpos(3)*0.04)^2; %marker width will be exactly 4% of 
mAlpha = 0.65;                       %the x-axis width.

%recovery plot
h = gobjects(size(parameters));
recoveryCounts  = [0 0];
for n = 1:length(parameters)
    %select a parameter
    p = parameters{n};
    %all parameter instances, determine X and Y coordinates
    trueVals = trueValues.(p)(:);
    try
        estVals = estimatedValues.(p).(pointEstType)(:);
    catch
        error('pointEstType ''%s'' was not present in estimatedValues.', ...
            pointEstType)
    end
    try
        lowerCIbound = estimatedValues.(p).(lowerBound)(:);
        upperCIbound = estimatedValues.(p).(upperBound)(:);
    catch
        error(['lower/upper bounds for credIntType ''%s'' were not ' ...
            'found in estimatedValues.'],credIntType)
    end
    % ... determine which points to plot are 
    inCI = lowerCIbound<=trueVals & trueVals<=upperCIbound;
    outsideCI = ~inCI;
    %plot 95% credible interval
    if any(outsideCI) && makeErrorBars
        g = plot([trueVals(outsideCI) trueVals(outsideCI)]', ...
            [lowerCIbound(outsideCI) upperCIbound(outsideCI)]', ...
            'marker','none','linewidth',lWidth,'color',failedColor);
        %transparent
        for ng = 1:length(g)
            g(ng).Color(4) = 0.5;
        end
    end
    if any(inCI) && makeErrorBars
        g = plot([trueVals(inCI) trueVals(inCI)]', ...
            [lowerCIbound(inCI) upperCIbound(inCI)]', ...
            'marker','none','linewidth',lWidth,'color',paramColors(n,:));
        %transparent
        for ng = 1:length(g)
            g(ng).Color(4) = 0.5;
        end
    end
    %plot marker for true vs. estimated value
    h(n) = scatter(trueVals,estVals,mPtSize, ...
        'markerfacecolor',paramColors(n,:),'markerfacealpha',mAlpha, ...
        'markeredgecolor','none');
    if any(outsideCI)
        %if failed to recover, add a red X
        scatter(trueVals(outsideCI),estVals(outsideCI),mPtSize*1.5, ...
            'marker','x','markeredgecolor',failedColor,'linewidth',lWidth*1.5);
%         %if failed to recover, add a red halo
%         scatter(trueVals(outsideCI),estVals(outsideCI),mPtSize*1.5, ...
%             'markerfacecolor','none','markeredgecolor',failedColor, ...
%             'linewidth',lWidth*2);
    end
    if exist('minBounds','var')
        xl = xlim; yl = ylim;
        xlim([min(minBounds(1),xl(1)) max(minBounds(2),xl(2))])
        ylim([min(minBounds(1),yl(1)) max(minBounds(2),yl(2))]) 
    end
    %keep a running count of how many are/are not in credible interval
    recoveryCounts = recoveryCounts + [sum(inCI==1) sum(inCI==0)];
end

%overlay a perfect recovery line
xl = ax.XLim; yl = ax.YLim;
recov(1) = min(xl(1),yl(1));
recov(2) = max(xl(2),yl(2));
plot(recov,recov,'k')

%optionally, overlay a line to highlight a meaningful value
if ~isempty(highlightValue)
    %redraw recovery if needed
    recov(1) = min(recov(1),highlightValue);
    recov(2) = max(recov(2),highlightValue);
    plot(recov,recov,'k:')
    %highlight value
    plot(highlightValue*[1 1],recov,'color','k','linewidth',lWidth,'linestyle',':')
    plot(recov,highlightValue*[1 1],'color','k','linewidth',lWidth,'linestyle',':')
end

% %add correlation text
% if length(parameters)==1
%     xl = xlim; yl = ylim;
%     text(xl(2)-+0.01*diff(xl),yl(1)+0.01*diff(yl),sprintf('rho = %.3f',corr(estVals,trueVals)), ...
%         'HorizontalAlignment','right','VerticalAlignment','bottom')
% end

%format recovery plot
% axis square
xlim(recov)
ylim(recov)
ticks = get(gca,'xtick'); ytick = get(gca,'ytick');
if ~isequal(ticks,ytick)
    if length(yticks) > length(ticks), ticks = yticks; end
end
set(gca,'xtick',ticks,'ytick',ticks)
%... make a static grid
ax = gca;
ax.XAxis.MinorTickValuesMode = 'manual';
ax.YAxis.MinorTickValuesMode = 'manual';
ax.XAxis.MinorTickValues = ticks;
ax.YAxis.MinorTickValues = ticks;
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
grid minor
%other formatting
xlabel('true value')
ylabel('estimated value')
set(gca,'fontsize',figFontSz)
if length(parameters) > 1
    legend(h,parameters,'interpreter','none', ...
        'location','eastoutside','orientation','vertical','box','off')
    fpos = f.Position;
    f.Position = [fpos(1:2) fpos(3)+100*figSzScaling fpos(4)];
else
    title(parameters{1},'interpreter','none')
end

%% output
if nargout == 1
    varargout{1} = recoveryCounts;
end

end
