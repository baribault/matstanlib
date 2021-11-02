function plotautocorr(varargin)
%PLOTAUTOCORR plots the estimated autocorrelation function for each chain.
% 
% this function generates a figure which, for a given parameter, presents
% the normalized autocorrelation function for the posterior samples in each
% chain as overlaid stem plots (with a small amount of jitter). 
% 
% 
% PLOTAUTOCORR(SAMPLES,PARAMETER)
%   SAMPLES is a structure containing posterior samples (consistent with
%   the format generated by extractsamples.m).  
%   PARAMETER is the name of a scalar parameter, or the (indexed) name 
%   of a parameter instance (e.g., 'mu[3]', 'sigma[2,1]').
%   for this syntax, SAMPLES and PARAMETER are both requrired inputs.
% 
% PLOTAUTOCORR(CHAINS)
%   alternatively, CHAINS is a [nIterations nChains]-sized matrix of
%   posterior samples, may be given.  
%   for this alternative synatx, CHAINS is the only required input.  
% 
% 
% PLOTAUTOCORR(...,MAXLAG)
%   in either syntax, MAXLAG may also be given, to determine the maximum
%   lag presented. the default value is 50.  
% 
% 
% See also PLOTLP, RHATTABLE, EXTRACTSAMPLES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

msl.options

%% parse required inputs
if nargin == 0
    error('at least one input is required.')
end

%%% SYNTAX #1 %%%
if isstruct(varargin{1})
    %first syntax option ...
    %   SAMPLES     = a structure of mcmc samples
    %   PARAMETER   = a parameter name string, which must be a valid 
    %                 scalar parameter name or parameter instance name
    numRequiredInputs = 2;
    if nargin < numRequiredInputs
        error(['if the first input is the samples structure, then ' ...
            'a second input (a parameter name string) must also be given.'])
    elseif ~ischar(varargin{2})
        error(['if the first input is the samples structure, then ' ...
            'the second input must be a parameter name string.'])
    end
    %parse parameter input
    [parameter,ind] = str2ind(varargin{2});
    if ~ismember(parameter,fieldnames(varargin{1}))
        error(['the parameter name ''%s'' was not found in the ' ...
            'samples structure.'],parameter)
    end
    if isempty(ind)
        %parameter name was given
        if ismatrix(varargin{1}.(parameter))
            %scalar parameter
            chains = varargin{1}.(parameter);
        else
            %not actually scalar parameter
            error(['please modify the parameter name input (%s) ' ...
                'to specify a parameter instance (e.g., %s[3]).'], ...
                parameter,parameter)
        end
    else
        %instance name was given
        instance = varargin{2};
        nParamDims = ndims(varargin{1}.(parameter)) - 2;
        if length(ind) ~= nParamDims
            error(['the number of indices (%i) in the given ' ...
                'parameter instance name %s does not match the ' ...
                'actual number of parameter dimensions (%i).'], ...
                length(ind),varargin{2},nParamDims)
        end
        chains = varargin{1}.(parameter)(:,:,ind{:});
        parameter = instance;
    end
%%% SYNTAX #2 %%%
elseif isnumeric(varargin{1})
    %second syntax option ...
    %   CHAINS     = a matrix of mcmc samples of size [nIterations nChains]
    numRequiredInputs = 1;
    if ismatrix(varargin{1})
        chains = varargin{1};
    else
        error('chains input must be of size [nIterations nChains].')
    end
    parameter = ''; %no parameter name
else
    error(['invalid number of required inputs or input type(s). ' ...
        'type ''help plotautocorr'' for info.'])
end

%% parse optional inputs
varargin = varargin(numRequiredInputs+1:end);

maxLag = 50;                foundMaxLag = false;
validPlotTypes = {'stem','bar','line'};
plotType = 'stem';          foundPlotType = false;

for v = 1:length(varargin)
    if isempty(varargin{v})
        %do nothing
    %plotType
    elseif ischar(varargin{v})
        if ~ismember(varargin{v},validPlotTypes)
            error(['''%s'' is not a valid autocorreltion plot type.  \n' ...
                'valid plot types include: ''%s'''], ...
                varargin{v},strjoin(validPlotTypes,''', '''))
        elseif foundPlotType
            error('cannot have two plotType inputs.')
        else
            plotType = varargin{v};
            foundPlotType = true;
        end
    %maxLag
    elseif isnumeric(varargin{v})
        if ~isscalar(varargin{v}) || mod(varargin{v},1)>0 || ~(varargin{v} > 1)
            error('maxLag must be a single integer > 1.')
        elseif varargin{v} > size(chains,1)-1
            error(['cannot compute the autocorrelation function ' ...
                'with the requested maximum lag (%i). the maximum ' ...
                'lag must be less than the number of iterations ' ...
                '(must be < %i).'],maxLag,size(chains,1))
        elseif foundMaxLag
            error('only one maxLag input may be given.')
        else
            maxLag = varargin{v};
            foundMaxLag = true;
        end
    else
        error('optional input type not recognized.')
    end
end

%% compute **normalized** autocorrelation function, chain by chain
nIterations = size(chains,1);
nChains = size(chains,2);

%adjust maxLag based on the number of iterations
maxLag = min(maxLag,nIterations-1);
nLags = min(maxLag+1,nIterations);

norm_acf = NaN([nLags nChains]);
for m = 1:nChains
    [auco,~] = acf(chains(:,m));
    norm_acf(:,m) = auco(1:nLags)/auco(1); %**normalized** autocorrelation
end

%% plot autocorrelation function
%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos(1:2) [600 400]*figScaling];
close(dumf.Number); %close dummy figure
hold on;

%get chain colors (matches tracedensity)
chainColors = getcolors('mcmc');

%autocorrelation plot
%underlay a line at rho == 0
line([-0.5 maxLag+0.5],[0 0],'color',0.8*[1 1 1],'linestyle',':')
%empirical acf
h = gobjects([nChains 1]);
switch plotType
    case 'stem'
        jitter = 0.2/nChains;
        x = (0:maxLag) - jitter*nChains/2;
        for m = 1:nChains
            h(m) = stem(x,norm_acf(:,m), ...
                'color',getcolors('lightgray'), ...
                'markerSize',round(markSz/2), ...
                'markerfacecolor',chainColors(m,:),'markeredgecolor',chainColors(m,:));
            x = x + jitter;
        end
    case 'line'
        x = 0:maxLag;
        for m = 1:nChains
            h(m) = plot(x,norm_acf(:,m), ...
                'color',chainColors(m,:),'linewidth',linePt);
        end
end
xlim([-0.5 maxLag+0.5])
xlabel('lag')
ylabel('autocorrelation function')
% set(gca,'ylim',[min(0,min(acf(:))-(1-min(acf(:)))*0.1) 1])
if ~isempty(parameter)
    title(parameter,'interpreter','none')
end
set(gca,'box','on','fontsize',fontSz)

%legend with chain numbers
legend(h,arrayfun(@num2str,1:nChains,'uni',0),'box','off')

end