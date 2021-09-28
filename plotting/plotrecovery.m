function varargout = plotrecovery(samples,trueValues,parameterRequest, ...
    varargin)
%PLOTRECOVERY plots true parameter values against estimated values.
%   
% this function visualizes the quality of parameter recovery.  
% when a Bayesian model is applied to simulated data, "recovery" is
% operationalized as whether or not each of the known parameter values used
% to simulate the data (called the "true" parameter values) are within the
% corresponding estimated credible interval. 
% 
% in each generated figure, the true values are plotted against the model's
% estimated values with semi-transparent markers, and a perfect recovery
% line is overlaid in black.  parameters that were not successfully
% recovered are indicated with a red X. 
% 
% 
% PLOTRECOVERY(SAMPLES,TRUEVALUES)
% PLOTRECOVERY(ESTIMATEDVALUES,TRUEVALUES)
%   generates a recovery plot for all parameters based on TRUEVALUES, a
%   structure of the true parameter values used used to simulate the data,
%   and *either* SAMPLES, the strucuture of posterior samples, *or*
%   ESTIMATEDVALUES, a structure of estimated values.  
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
% optional inputs may be given as property name-value pairs:
% PLOTRECOVERY(SAMPLES,TRUEVALUES,'name1',VALUE1,'name1',VALUE2,...)
% PLOTRECOVERY(ESTIMATEDVALUES,TRUEVALUES,'name1',VALUE1,'name1',VALUE2,...)
% pairs may be given in any order.  property names are case insensitive.
%   
%   'pointesttype'
%       POINTESTTYPE is a string indicating what posterior statitistic
%       should be used as the point parameter estimate.  
%       valid options include: 
%           'mean'    >>  mean of the posterior samples
%           'median'  >>  median of the posterior samples
%           'mode'    >>  mode of the posterior samples [for discrete only]
%           'argmax'  >>  maximum a posteriori
%       the default value for POINTESTTYPE is 'mean'.
%   
%   'markercolor'
%       MARKERCOLOR can offer a [nParameters 3] matrix of RGB-01 colors,
%       one for each recovered parameter. 
%   
%   'singleFig'
%       if SINGLEFIG is true/1, then a single figure will be generated that
%       presents the recovery of *all* parameters, with a legend.
%       if SINGLEFIG is false/0, then one figure per parameter will be
%       generated.  each figure presents the recovery of all instances of a
%       given parameter, denoted in the title. 
%       the default value for SINGLEFIG is true.
%   
%   'CImethod'
%       CIMETHOD indicates how the credible interval (which is used to
%       determine whether the true value is or is not within the CI) should
%       be computed.  valid CIMETHOD strings include:
%           'central'   >>  central (i.e., equal-quantile-tailed) interval
%           'HDI'       >>  highest density interval (HDI) [continuous]
%       the default value is 'central'.
%   
%   'CIwidth'
%       CIWIDTH is a scalar numeric value between 0 and 1 indicating what
%       probability credible interval should used be used to determine
%       whether or not the true value is within the CI.
%       the default value for CIWIDTH is 0.95 (for 95% credible intervals).
%   
%   'addinterval'
%       ADDINTERVAL may be true or false.
%       if ADDINTERVAL is true, then credible intervals are displayed on
%       the plot as capless "error bars".  
%       if successfully recovered (the true value is in the CI), then the
%       CI color will match the marker color.   
%       if failed to recover (the true value is outside the CI), then the
%       CI color will be red.   
%       the default value for SHOWINTERVAL is false (error bars are omitted).
%   
%   'criticalvalue'
%       CRITICALVALUE is a single numeric value.  
%       if CRITICALVALUE is given, a vertical line & a horizontal line
%       will be overlaid at the requested value (such as at 0).
%       if CRITICALVALUE is empty, no value is highlighted. 
%       by default, no value is highlighted.
%   
%   'addcorrelation'
%       ADDCORRELATION triggers the calculation and overlay of the
%       text reporting correlation between the true and estimated values.
%       *** NOTE: the correlation is a simple call to corr.m --- it is NOT
%                 a Bayesian correlation!
% 
% 
% See also GETSAMPLESTATS, COLLECTTRUEVALUES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

msl.options

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
if nargin < 2
    error('at least two inputs are required.')
end

%samples/estimatedValues
if isstruct(samples)
    inputFields = fieldnames(samples); 
    inputSz = size(samples.(inputFields{1}),[1 2]); %[nIters nChains]
    if isnumeric(samples.(inputFields{1}))
      if all(structfun(@(x) isequal(size(x,[1 2]),inputSz),samples))
        %yep, it's samples
        gaveSamples = true;        
      else
        %nope, it's trueValues
        error(['first input must be either the samples structure ' ...
            '(in the format generated by extractsamples.m) ' ...
            'or an estimatedValues structure ' ...
            '(in the format generated by getsamplestats.m).  ' ...
            'did you input the trueValues structure first by mistake?'])
      end
    elseif isstruct(samples.(inputFields{1})) ...
            && isfield(samples.(inputFields{1}),'mean')
        %it's actually estimatedValues
        gaveSamples = false; 
        estimatedValues = samples; clearvars samples
    else
        error(['first input must be either the samples structure ' ...
            '(in the format generated by extractsamples.m) ' ...
            'or an estimatedValues structure ' ...
            '(in the format generated by getsamplestats.m).'])
    end
else
    error(['first input must be either the samples structure ' ...
        '(in the format generated by extractsamples.m) ' ...
        'or an estimatedValues structure ' ...
        '(in the format generated by getsamplestats.m).'])
end

%trueValues
if ~isstruct(trueValues)
    error('the second input must be a structure of true parameter values.')
end

%% determine which subset of parameters to plot
%get parameter names
if gaveSamples
    allKnownParameters = union(fieldnames(samples),fieldnames(trueValues));
    matchedParameters = intersect(fieldnames(samples),fieldnames(trueValues));
else
    allKnownParameters = union(fieldnames(estimatedValues),fieldnames(trueValues));
    matchedParameters = intersect(fieldnames(estimatedValues),fieldnames(trueValues));
end

%parse parameterRequest
if nargin < 3 || isempty(parameterRequest)
    parameterRequest = matchedParameters; %default
elseif ischar(parameterRequest) 
    if ismember(parameterRequest,allKnownParameters)
        parameterRequest = {parameterRequest};
    else
        %actually just a property name!!!
        varargin = [{parameterRequest},varargin];
        parameterRequest = matchedParameters; %default
    end
elseif iscell(parameterRequest) && all(cellfun(@ischar,parameterRequest))
    %all good
else
    error(['the third input must be parameterRequest.  the requested ' ...
        'parameter name(s) must be input as a string or cell of strings ' ...
        '(or an empty cell, to request all parameters).'])
end

%check all requested parameter names are valid
if any(~ismember(parameterRequest,matchedParameters))
    missingRequests = setdiff(parameterRequest,matchedParameters);
    if gaveSamples, str = 'samples'; else, str = 'estimatedValues'; end
    warning(['the following elements of parameterRequest were ' ...
        'not fields in the %s structure and/or ' ...
        'were not fields in the trueValues structure: %s'], ...
        str,['''' strjoin(missingRequests,''', ''') ''''])
    parameterRequest = setdiff(parameterRequest,missingRequests);
end

%make sure there's still some parameters left!!!
if isempty(parameterRequest)
    error('no valid parameters available to plot!')
end

% ... just rename at this point
parameters = parameterRequest;
nParameters = length(parameters);

%% parse optional inputs
%ensure all remaining inputs are consistent with name-value pair format
if mod(length(varargin),2) > 0 || ...
        ~all(cellfun(@ischar,varargin(1:2:length(varargin))))
    error(['all inputs (after samples and parameterRequest) ' ...
        'must be submitted as name-value pairs.'])
elseif length(unique(varargin(1:2:length(varargin))))<length(varargin(1:2:length(varargin)))
    error('at least one property name (in the name-value pair inputs) is a duplicate.')
end

%set default values for the optional inputs
pointEstType = 'mean';      validPointEstTypes = {'mean','median','mode','argmax'};
markerColor = [];       %use default colors
singleFig = true;       %default
CImethod = 'central';       validCIMethods = {'central','hdi'};
CIwidth = 0.95;
addInterval = false;    %default
criticalValue = [];     %default
addCorrelation = false; %default

varargin(1:2:length(varargin)) = ... %case insensitive (convert all to lower case)
    cellfun(@lower,varargin(1:2:length(varargin)),'uni',0);
for v = 1:2:length(varargin)
    switch varargin{v}
        %----------------------------------------------------------------%
        case 'pointesttype'
            pointEstType = lower(varargin{v+1});
            if ~ischar(pointEstType)
                error('pointEstType must be a string (@ischar==true).')
            elseif ~ismember(pointEstType,validPointEstTypes)
                error(['pointEstType input is not valid. valid inputs ' ...
                    'include: ''%s''.'],strjoin(validPointEstTypes,''', '''))
            end
        %----------------------------------------------------------------%
        case 'markercolor'
            markerColor = varargin{v+1};
            if isempty(markerColor)
                %do nothing --- will use default colors
            elseif ischar(markerColor)
                markerColor = getcolors(lower(markerColor));
            elseif isnumeric(markerColor) 
                if ~any(isnan(markerColor(:))) || ...
                        ~(all(markerColor(:) >= 0) && all(markerColor(:) <= 1))
                    error(['markerColor must be a valid RGB-01 color ' ...
                        'but values were given outside the [0 1] range ' ...
                        'or were NaN.'])
                elseif ismatrix(markerColor) && size(markerColor,2)==3
                    %one or more colors
                    if size(markerColor,1) < nParameters
                        error(['the number or rows in markerColor is less than ' ...
                            'the number of parameters.'])
                    end
                else
                    error(['markerColor must be a Nx3 vector or matrix, ' ...
                        'but the size given was %s.'],mat2str(size(markerColor)))
                end
            else
                error('markerColor must be an RGB-01 color or a known color name.')
            end
        %----------------------------------------------------------------%
        case 'singlefig'
            singleFig = varargin{v+1};
            if ~isscalar(singleFig) || ~ismember(singleFig,[0 1])
                error('singleFig input must be true, 1, false, or 0.')
            end
        %----------------------------------------------------------------%
        case 'cimethod'
            CImethod = lower(varargin{v+1});
            if ~ischar(CImethod)
                error('CImethod must be a string (@ischar==true).')
            elseif ~ismember(CImethod,validCIMethods)
                error(['CImethod input is not valid. valid inputs ' ...
                    'include: ''%s''.'],strjoin(validCImethods,''', '''))
            end
        %----------------------------------------------------------------%
        case 'ciwidth'
            CIwidth = varargin{v+1};
            if ~isnumeric(CIwidth) || ~isscalar(CIwidth)
                error('CIwidth must be a single number.')
            elseif CIwidth < 0 || CIwidth > 1
                error('CIwidth must be a proportion [0,1].')
            end
        %----------------------------------------------------------------%
        case 'addinterval'
            addInterval = varargin{v+1};
            if ~isscalar(addInterval) || ~ismember(addInterval,[0 1])
                error('addInterval input must be true, 1, false, or 0.')
            end
        %----------------------------------------------------------------%
        case 'criticalvalue'
            criticalValue = varargin{v+1};
            if ~isnumeric(criticalValue) || ~isscalar(criticalValue)
                error('criticalvalue must be a single number.')
            end
        %----------------------------------------------------------------%
        case 'addcorrelation'
            addCorrelation = varargin{v+1};
            if ~isscalar(addCorrelation) || ~ismember(addCorrelation,[0 1])
                error('addCorrelation input must be true, 1, false, or 0.')
            end
        %----------------------------------------------------------------%
        otherwise
            error('the optional input name ''%s'' was not recognized.', ...
                varargin{v})
    end
end

%% prepare to plot
%get some colors
recoveredColor = getcolors('gray');    %cred int bars
failedColor = getcolors('red');        %cred int bars
if ~isempty(markerColor)
    if isvector(markerColor)
        %same color for all parameters
        paramColors = repmat(markerColor,[nParameters 1]);
    else
        %different color for each parameter
        paramColors = markerColor;
    end
else
    paramColors = getcolors('ocean');           %point estimate markers
    if nParameters >= size(paramColors,1)
        paramColors = getcolors('mcmc');
        if nParameters >= size(paramColors,1)
            %so many parameters omg why
            paramColors = rand([nParameters 3]);
        end
    end
end

%and pick your CI type
if gaveSamples
    %do nothing
else
    switch CImethod
        case 'central'
            CIstr = 'CI';
        case {'hdi','hpd','hpdi'}
            CIstr = 'HDI';
    end
    if ~ismember(CIwidth,[0.95,0.90,0.50])
        error(['if estimatedValues is given, the CIwidth can only be ' ...
            '0.95, 0.90, or 0.50.'])
    end
end

%% recovery plot
recoveryCounts  = [0 0];

%recovery plot
h = gobjects(size(parameters));
for p = 1:nParameters
    %one figure or multiple figures?
    if (singleFig && p==1) || ~singleFig
        %start a figure
        dumf = figure(999); %dummy figure to protect sizing
        f = figure('color',[1 1 1]);
        fpos = f.Position;
        f.Position = [fpos([1 2]) [425 425]*figScaling];
        close(dumf.Number); %close dummy figure
        %start an axis
        ax = axes;
        hold on
        
        %set marker size in pts (just the once)
        if p==1
            %marker size ... sigh
            %temporarily mess with plot
            axunits = get(gca,'units'); %original settings
            set(gca,'units','points');  %temporarily change settings
            axpos = get(gca,'pos');
            set(gca,'units',axunits);   %revert to original settings
            %plot settings
            mPtSize = (axpos(3)*0.04)^2;%marker width will be exactly 4% of 
            mAlpha = 0.8;%0.65;         %the x-axis width.
        end
    end
    
    %--------------------------------------------------------------------%
    %select a parameter
    parameter = parameters{p};
    %all parameter instances, determine X coordinates: true values
    trueVals = trueValues.(parameter)(:);
    %all parameter instances, determine Y coordinates: estimated values
    if gaveSamples
        %%% from samples %%%
        chainSize = size(samples.(parameter));
        if length(chainSize) == 2
            %scalar parameter
            parameterSize = [1 1];
        elseif length(chainSize) == 3
            %vector parameter
            parameterSize = [chainSize(end) 1];
        else
            %N-dimensional matrix parameter
            parameterSize = chainSize(3:end);
        end
        instanceNames = getparaminstances(parameter, ...
            fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));
        nInstances = length(instanceNames);
        %compute point estimate & credible interval
        estVals = NaN(parameterSize);
        lowerCIbound = NaN(parameterSize);
        upperCIbound = NaN(parameterSize);
        for n = 1:nInstances
            [ind{1:length(parameterSize)}] = ind2sub(parameterSize,n);
            if nInstances == 1
                chains = samples.(parameter);
            else
                chains = samples.(parameter)(:,:,ind{:});
            end
            chains = chains(:);
            switch pointEstType
                case 'mean',    estVals(ind{:}) = mean(chains);
                case 'median',  estVals(ind{:}) = median(chains);
                case 'mode',    estVals(ind{:}) = mode(chains);
                case 'argmax' 
                    [f,x] = smoothdensity(chains);
                    estVals(ind{:}) = x(f==max(f));
            end
            CI = computecredint(chains,CImethod,CIwidth);
            lowerCIbound(ind{:}) = CI(1);
            upperCIbound(ind{:}) = CI(2);
        end
    else
        %%% from estimatedValues %%%
        try
            estVals = estimatedValues.(parameter).(pointEstType)(:);
        catch
            error('pointEstType ''%s'' was not present in estimatedValues.', ...
                pointEstType)
        end
        try
            lowerCIfield = ['lower' CIstr num2str(CIwidth*100)];
            lowerCIbound = estimatedValues.(parameter).(lowerCIfield);
            upperCIfield = ['upper' CIstr num2str(CIwidth*100)];
            upperCIbound = estimatedValues.(parameter).(upperCIfield);
        catch
            error(['the fields ''lower%s'' and/or ''upper%s'' were not ' ...
                'found in estimatedValues.%s.'],lowerCIfield,upperCIfield,parameter)
        end
    end
    % ... determine which points to plot are 
    inCI = lowerCIbound<=trueVals & trueVals<=upperCIbound;
    outsideCI = ~inCI;
    %plot 95% credible interval
    if any(outsideCI) && addInterval
        g = plot([trueVals(outsideCI) trueVals(outsideCI)]', ...
            [lowerCIbound(outsideCI) upperCIbound(outsideCI)]', ...
            'marker','none','color',failedColor);
            % 'marker','none','linewidth',linePt,'color',failedColor);
        %transparent
        for ng = 1:length(g)
            g(ng).Color(4) = 0.5;
        end
    end
    if any(inCI) && addInterval
        g = plot([trueVals(inCI) trueVals(inCI)]', ...
            [lowerCIbound(inCI) upperCIbound(inCI)]', ...
            'marker','none','linewidth',linePt,'color',paramColors(p,:));
        %transparent
        for ng = 1:length(g)
            g(ng).Color(4) = 0.5;
        end
    end
    %plot marker for true vs. estimated value
    h(p) = scatter(trueVals,estVals,mPtSize, ...
        'markerfacecolor',paramColors(p,:),'markerfacealpha',mAlpha, ...
        'markeredgecolor','none');
    if any(outsideCI)
        %if failed to recover, add a red X
        scatter(trueVals(outsideCI),estVals(outsideCI),mPtSize*1.5, ...
            'marker','x','markeredgecolor',failedColor,'linewidth',linePt*1.5);
%         %if failed to recover, add a red halo
%         scatter(trueVals(outsideCI),estVals(outsideCI),mPtSize*1.5, ...
%             'markerfacecolor','none','markeredgecolor',failedColor, ...
%             'linewidth',linePt*2);
    end
    if exist('minBounds','var')
        xl = xlim; yl = ylim;
        xlim([min(minBounds(1),xl(1)) max(minBounds(2),xl(2))])
        ylim([min(minBounds(1),yl(1)) max(minBounds(2),yl(2))]) 
    end
    %keep a running count of how many are/are not in credible interval
    recoveryCounts = recoveryCounts + [sum(inCI==1) sum(inCI==0)];
    %--------------------------------------------------------------------%
    
    %one figure or multiple figures?
    if (singleFig && p==nParameters) || ~singleFig
        %overlay a perfect recovery line
        xl = ax.XLim; yl = ax.YLim;
        recov(1) = min(xl(1),yl(1));
        recov(2) = max(xl(2),yl(2));
        plot(recov,recov,'k')
        
        %optionally, overlay a line to highlight a meaningful value
        if ~isempty(criticalValue)
            %redraw recovery if needed
            recov(1) = min(recov(1),criticalValue);
            recov(2) = max(recov(2),criticalValue);
            plot(recov,recov,'k:')
            %highlight value
            plot(criticalValue*[1 1],recov,'color','k','linewidth',linePt,'linestyle',':')
            plot(recov,criticalValue*[1 1],'color','k','linewidth',linePt,'linestyle',':')
        end
        
        %add correlation text
        if addCorrelation && ((singleFig && nParameters==1) || ~singleFig)
            xl = xlim; yl = ylim;
            % text(xl(1)+0.5*diff(xl),yl(2)-0.01*diff(yl), ...
            text(xl(2)-+0.01*diff(xl),yl(1)+0.01*diff(yl), ...
                sprintf('\\rho = %.3f',corr(estVals,trueVals)), ...
                'interpreter','tex','fontsize',fontSz, ...
                'HorizontalAlignment','right','VerticalAlignment','bottom')
                % 'HorizontalAlignment','center','VerticalAlignment','top')
        end
        
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
        %... other formatting
        xlabel('true value')
        ylabel('estimated value')
        set(gca,'fontsize',fontSz)
        if singleFig && nParameters > 1
            legend(h,parameters,'interpreter','none', ...
                'location','eastoutside','orientation','vertical','box','off')
            fpos = f.Position;
            f.Position = [fpos(1:2) fpos(3)+100*figScaling fpos(4)];
        else
            title(parameter,'interpreter','none')
        end
    end
end

%% output
if nargout == 1
    varargout{1} = recoveryCounts;
end

end
