function jointdensity(samples,parameter1,parameter2,varargin)
%JOINTDENSITY plots a marginal bivariate posterior density.
% 
% this function generates a single figure with three panels:
%   (1) a scatter plot of the samples for PARAMETER1 versus the samples 
%       for PARAMETER2, 
%   (2) a smoothed univariate marginal density for PARAMETER1, and
%   (3) a smoothed univariate marginal density for PARAMETER2.
% 
% 
% JOINTDENSITY(SAMPLES,PARAMETER1,PARAMETER2)
%   SAMPLES is a structure of posterior samples (in the format generated by
%   extractsamples.m).  
%   PARAMETER1 and PARAMETER2 are names of the scalar parameters or
%   parameter instance names to be plotted. 
% 
% 
% JOINTDENSITY(...,DIAGNOSTICS)
% JOINTDENSITY(...,DIVERGENT)
%   optionally, samples derived from divergent transitions may be 
%   highlighted in red on the plot.  
%   either submit a DIAGNOSTICS structure directly, or submit a 
%   [nIterations nChains]-sized matrix of indicators of divergent 
%   transitions.
% 
% 
% Examples:
%   jointdensity(samples,'mu_group','sigma[1]',diagnostics)
%   jointdensity(samples.mu_group,samples.sigma(:,:,1), ...
%       diagnostics.divergent__)
%   
% See also SMOOTHDENSITY, PLOTDENSITY, MULTIDENSITY, VIOLINDENSITY, 
% 	EXTRACTSAMPLES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

matstanlib_options

%% parse inputs
if nargin < 3
    error('at least three inputs are required.')
end

%samples
if ~isstruct(samples)
    error('first input must be the samples structure.')
end

%parameter1, parameter2
if ~ischar(parameter1), error('parameter1 must be a string (@ischar==true).'); end
if ~ischar(parameter2), error('parameter2 must be a string (@ischar==true).'); end
for n = 1:2
    %... don't want to repeat code ...
    if n==1,     [parameter,ind] = str2ind(parameter1); paramInput = parameter1;
    elseif n==2, [parameter,ind] = str2ind(parameter2); paramInput = parameter2;
    end
    if isempty(ind)
        %parameter name given
        if ismatrix(samples.(parameter))
            %scalar parameter
            if n==1,     chains1 = samples.(parameter)(:);
            elseif n==2, chains2 = samples.(parameter)(:);
            end
        else
            %not actually scalar parameter
            error(['please modify the parameter%i input (''%s'') ' ...
                'to specify which parameter instance (e.g., ''%s[3]'') ' ...
                'is to be used.'],n,parameter,parameter)
        end
    else
        %instance name given
        nParamDims = ndims(samples.(parameter)) - 2;
        if length(ind) ~= nParamDims
            error(['the number of indices (%i) in the given parameter%i ' ...
                'instance name does not match the actual number of ' ...
                'parameter dimensions (%i).'],length(ind),n,nParamDims)
        end
        if n==1,     chains1 = reshape(samples.(parameter)(:,:,ind{:}),[],1);
        elseif n==2, chains2 = reshape(samples.(parameter)(:,:,ind{:}),[],1);
        end
    end
end

%%% optional input %%%
if isempty(varargin)
    %do nothing
elseif length(varargin)==1
    %highlight divergent transitions
    if isstruct(varargin{1}) && isfield(varargin{1},'divergent__')
        divergent = varargin{1}.divergent__;
    elseif isnumeric(varargin{1}) && ismatrix(varargin{1})
        divergent = varargin{1};
    else
        error(['the optional input must be either the diagnostics ' ...
            'structure or a matrix of divergent transition indicators.'])
    end
else
    error('too many inputs.')
end

% chains1 = binornd(10,0.3,size(chains1)); %test discrete

%% plot bivariate density
%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = get(f,'pos');
f.Position = [fpos(1:2) [575 575]*figScaling];
close(dumf.Number) %close dummy figure
%... and three axes
% axJ  = axes('pos',[0.1 0.1 0.6 0.6]);
axJ  = axes('pos',[0.1 0.1 0.58 0.58]);
axM1 = axes('pos',[0.1 0.75 0.58 0.2]);
axM2 = axes('pos',[0.75 0.1 0.2 0.58]);

%... and a nice color
sampleColor = getcolors('gray');
sampleAlpha = 0.125;
divergentColor = getcolors('red');
divergentAlpha = 0.65;

%(1) plot joint density
axes(axJ); hold on
scatter(chains1,chains2,'sizedata',markSz, ...
    'markerfacecolor',sampleColor,'markeredgecolor',sampleColor, ...
    'markerfacealpha',sampleAlpha,'markeredgealpha',sampleAlpha)
if exist('divergent','var')
    scatter(chains1(divergent==1),chains2(divergent==1),'sizedata',markSz, ...
    'markerfacecolor',divergentColor,'markeredgecolor',divergentColor, ...
    'markerfacealpha',divergentAlpha,'markeredgealpha',divergentAlpha)
end
%format plot
set(gca,'fontsize',fontSz,'box','on','xgrid','on','ygrid','on')
if exist('parameter1','var'); xlabel(parameter1,'fontweight','bold','interpreter','none'); end
if exist('parameter2','var'); ylabel(parameter2,'fontweight','bold','interpreter','none'); end
% axis square %NO! this throws the alignment with other axes out of whack

%(2) plot marginal #1
axes(axM1)
plotmarginal(axM1,chains1,sampleColor,linePt)
if exist('parameter1','var'); title(parameter1,'interpreter','none'); end
set(gca,'fontsize',fontSz)

%(3) plot marginal #2
axes(axM2)
plotmarginal(axM2,chains2,sampleColor,linePt)
if exist('parameter2','var'); title(parameter2,'interpreter','none'); end
set(gca,'fontsize',fontSz)
view([90 -90])

%check whether joint density scatter & marginals axis limits & ticks match
matchlimitsandticks(axJ,axM1,axM2)

end

%-----------------------------------------------------------------------%
function plotmarginal(ax,chains,color,linePt)
fcolor = (color + 2*[1 1 1])/3;
ecolor = color;
isDiscrete = all(~mod(chains,1));
[f,x] = smoothdensity(chains);
if isDiscrete
    %discrete plot
    bar(x,f,'facecolor',fcolor,'edgecolor',ecolor,'linewidth',linePt)
else
    %continuous plot
    area(ax,x,f,'facecolor',fcolor,'edgecolor',ecolor,'linewidth',linePt)
end
%format, either type of plot
set(ax,'box','off','xgrid','on')
end

%-----------------------------------------------------------------------%
function matchlimitsandticks(axJ,axM1,axM2)
%joint Y-axis and marginal #1 X-axis
%limits
if ~isequal(axJ.XLim,axM1.XLim)
    Jlim = axJ.XLim;
    Mlim = axM1.XLim;
    newlim = [min(Jlim(1),Mlim(1)) max(Jlim(2),Mlim(2))];
    axJ.XLim = newlim;
    axM1.XLim = newlim;
end
%ticks
if ~isequal(axJ.XTick,axM1.XTick)
    Jtick = axJ.XTick;
    Mtick = axM1.XTick;
    if length(Jtick) > length(Mtick)
        newtick = Jtick;
    elseif length(Jtick) < length(Mtick)
        newtick = Mtick;
    end
    axJ.XTick = newtick;
    axM1.XTick = newtick;
end
%joint Y-axis and marginal #2 X-axis
%limits
if ~isequal(axJ.YLim,axM2.XLim)
    Jlim = axJ.YLim;
    Mlim = axM2.XLim;
    newlim = [min(Jlim(1),Mlim(1)) max(Jlim(2),Mlim(2))];
    axJ.YLim = newlim;
    axM2.XLim = newlim;
end
%ticks
if ~isequal(axJ.YTick,axM2.XTick)
    Jtick = axJ.YTick;
    Mtick = axM2.XTick;
    if length(Jtick) > length(Mtick)
        newtick = Jtick;
    elseif length(Jtick) < length(Mtick)
        newtick = Mtick;
    end
    axJ.YTick = newtick;
    axM2.XTick = newtick;
end

end
