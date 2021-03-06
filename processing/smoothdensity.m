function [f,x] = smoothdensity(varargin)
%SMOOTHDENSITY uses posterior samples to approximate a smoothed density.
%  
% the function will first attempt to determine if the parameter is 
% discrete (i.e., integer-valued) or continuous.  then, the function 
% will approximate the posterior density based on the frequency of 
% posterior samples.  
% 
% NOTE: ideally, this function will call ksdensity.m, but if it is not 
% available (i.e, if MATLAB's Statistics Toolbox is not installed) a simple
% moving average will be used for smoothing.  
% 
% 
% [F,X] = SMOOTHDENSITY(SAMPLES,PARAMETER)
%   SAMPLES is a structure of posterior samples (in the format generated by
%   extractsamples.m). 
%   PARAMETER is the name of a scalar parameter or an indexed instance of a
%   nonscalar parameter (e.g., 'mu[3]','sigma[2,1]'). 
%   
%   kernel smoothing will be used to approximate the posterior
%   marginal density, F, over values of PARAMETER, X.
%  
% [F,X] = SMOOTHDENSITY(CHAINS)
%   alternatively, a single input, CHAINS, may be given.
%   CHAINS is a [nIterations nChains]-sized matrix of posterior samples.
%   this syntax is especially useful when calling the function from within
%   another function, to avoid storing a copy of SAMPLES in memory. 
% 
% 
% [F,X] = SMOOTHDENSITY(SAMPLES,PARAMETER,[NAME],[VALUE],...)
% [F,X] = SMOOTHDENSITY(CHAINS,[NAME],[VALUE],...)
%   in either syntax, options may also be given as name-value pairs, 
%   including: 'kernel','bandwidth','support'.  these are only used if
%   ksdensity.m is available; otherwise they are ignored.  see the help
%   file for ksdensity for details on these options and their defaults.
%   (note that the default for 'support' is based on the values in CHAINS.)
% 
% Examples, for a scalar parameter:
%   [F,X] = SMOOTHDENSITY(samples,'delta')
%   [F,X] = SMOOTHDENSITY(samples.delta)
% 
% Examples, for a non-scalar parameter:
%   [F,X] = SMOOTHDENSITY(samples,'sigma[2]')
%   [F,X] = SMOOTHDENSITY(samples.sigma(:,:,2))
% 
% 
% See also KSDENSITY, PLOTDENSITY
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% parse required inputs
if nargin == 0
    error('at least one input is required.')
end
if isstruct(varargin{1})
    numRequiredInputs = 2;
    if nargin < numRequiredInputs
        error('too few required inputs for this syntax.')
    end
    if ~ischar(varargin{2})
        error('at least one required input is an invalid type.')
    end
    %first syntax option ...
    %   SAMPLES     = a structure of mcmc samples
    %   PARAMETER   = a parameter name string, which must be a valid 
    %                 scalar parameter name or parameter instance name
    [parameter,ind] = str2ind(varargin{2});
    if isempty(ind)
        %parameter name given
        parameter = varargin{2};
        if ismatrix(varargin{1}.(parameter))
            %scalar parameter
            chains = varargin{1}.(parameter)(:);
        else
            %not actually scalar parameter
            error(['please modify the parameter name input (%s) ' ...
                'to specify which parameter instance (e.g., %s[3])' ...
                'is to be smoothed.'],parameter,parameter)
        end
    else
        %instance name given
        nParamDims = ndims(varargin{1}.(parameter)) - 2;
        if length(ind) ~= nParamDims
            error(['the number of indices (%i) in the given ' ...
                'parameter instance name %s does not match the ' ...
                'actual number of parameter dimensions (%i).'], ...
                length(ind),varargin{2},nParamDims)
        end
        chains = varargin{1}.(parameter)(:,:,ind{:});
        chains = chains(:);
    end
elseif isnumeric(varargin{1})
    numRequiredInputs = 1;
    %second syntax option ...
    %   CHAINS     = a matrix of mcmc samples of size [nIterations nChains]
    if ismatrix(varargin{1})
        chains = varargin{1}(:);
    else
        error('chains input must be of size [nIterations nChains].')
    end
else
    error(['invalid number of inputs or input type(s). ' ...
        'type ''help smoothdensity'' for info.'])
end

%% parse optional inputs
if nargin > numRequiredInputs
    validOptions = {'kernel','bandwidth','support'};
    options = varargin(numRequiredInputs+1:end);
    if mod(length(options),2) == 0
        %check if names in name-value pairs are valid options
        if any(~ismember(options(1:2:length(options)),validOptions))
            error('invalid optional input.')
        end
    else
        error('invalid number of optional inputs.')
    end
end

%% approximate a smooth posterior density based on the mcmc samples
isDiscrete = all(mod(chains,1)==0);
if isDiscrete
    %discrete parameter
    [f,~] = histcounts(chains,min(chains)-0.5:1:max(chains)+0.5);
    x = min(chains):max(chains);                %midpoints of bins
    f = f/sum(f);                               %normalize counts for pdf
else
    %continuous parameter
    if exist('ksdensity','file')
        %guess at the support
        mins = min(chains);
        maxs = max(chains);
        if 0 < mins && maxs < 1
            support = [0 1];
        elseif 0 < mins
            support = 'positive';
        else
            support = 'unbounded';
        end
        %ksdensity (requires stats toolbox)
        if exist('options','var')
            if ismember('support',options(1:2:length(options)))
                [f,x] = ksdensity(chains,options{:});
            else
                [f,x] = ksdensity(chains,'support',support,options{:});
            end
        else
            [f,x] = ksdensity(chains,'support',support);
        end
    else
        warning(['ksdensity.m not found. performing a simpler '
                'approximation with a moving average filter ...'])
        %well then i'll just do it myself!
        nbins = max(round(length(chains)/100/2),50);
        [f,edges] = histcounts(chains,nbins);
        binwidth = diff(edges(1:2));
        x = (edges(1:end-1) + edges(2:end))/2;  %midpoints of bins
        f = f/(sum(f)*binwidth);                %normalize counts for pdf
        oldnorm = sum(f);
        %smooth via moving average
        n = round(nbins/10);
        tmpy = [zeros([1 n]) f zeros([1 n])];
        for c = 1:length(f)
            f(c) = mean(tmpy(c:c+n*2));
        end
        f = f/sum(f)*oldnorm;
    end
end

end
