function MCSE = computemcse(chains,method,alpha)
%COMPUTEMCSE computes the Monte Carlo standard error (MCSE) from posterior samples. 
% 
% MCSE = COMPUTEMCSE(CHAINS)
% MCSE = COMPUTEMCSE(CHAINS,TYPE)
% MCSE = COMPUTEMCSE(CHAINS,'quantile',ALPHA)
%   this function computes the estimated Monte Carlo standard error (MCSE)
%   from CHAINS, a [nIterations nChains]-sized matrix of posterior samples,
%   according to current best computational practice.
%   
%   TYPE determines the quantity for which a MCSE should be computed.
%   many current MCSE types are available:
%       'mean'      >>  MCSE for the mean
%       'sd'        >>  MCSE for the standard deviation
%       'median'    >>  MCSE for the median
%     * 'quantile'  >>  MCSE for a particular quantile, ALPHA
% 
%   if TYPE is 'quantile', then an additional input, ALPHA, is required.  
%   ALPHA is the quantile for which the MCSE should be computed, given as a
%   proportion. 
%   as such, ALPHA should be a single numeric value between 0 and 1.
% 
% 
% Reference:  Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020, May). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% 
% 
% See also PLOTMCSE
% 
% (c) beth baribault 2021 ---                                 > matstanlib


%% check inputs
%chains
if nargin < 1, error('a chains matrix is required.'); end
if ~isnumeric(chains) && ~ismatrix(chains)
    error(['the first input must be chains, a [nIterations nChains]-sized ' ...
        'matrix of posterior samples.'])
end

%method
validMethods = {'mean','median','sd','quantile'};
if nargin < 2
    %type is REQUIRED!
    error(['a MCSE type or computation method is required. ' ...
        'valid computation methods include: ''%s'''], ...
        strjoin(validMethods,''', '''))
elseif ischar(method)
    if ~ismember(method,validMethods)
        error(['''%s'' is not a valid MCSE computation method.  \n' ...
            'valid computation methods include: ''%s'''], ...
            strjoin(validMethods,''', '''))
    end
else
    error('MCSE computation method must be a string (@ischar==true).')
end

%alpha
if isequal(method,'quantile')
    if nargin < 3
        error(['if MCSE type is ''quantile'', then alpha (the ' ...
            'requested MCSE quantile, as a proportion), is required.']) 
    elseif ~isnumeric(alpha) || ~isscalar(alpha)
        error('alpha must be a single numeric value.')
    elseif alpha <= 0 || alpha >= 1
        error('alpha must be in the interval (0,1).')
    end
end


%% compute MCSE

switch method
    % ------------------------------------------------------------------ %
    case 'mean'
        ESS_mean = computeess(chains,'mean');
        SD = std(chains(:),0); %normalize by N-1
        MCSE = SD / sqrt(ESS_mean);
    case 'sd'
        ESS_sd = computeess(chains,'sd');
        SD = std(chains(:),0); %normalize by N-1
        MCSE = SD * sqrt(exp(1)*((1 - 1/ESS_sd)^(ESS_sd-1)) - 1);
    case {'quantile','median'}
        if strcmp(method,'median'), alpha = 0.5; end
        ESS_quantile = computeess(chains,'quantile',alpha);        
        ab = betainv([0.1586553 0.8413447], ...
            ESS_quantile*alpha + 1,ESS_quantile*(1-alpha) + 1);
        sortedSamples = sort(chains(:));
        S = length(sortedSamples);
        thetaA = sortedSamples(max([floor(S*ab(1)) 1]));
        thetaB = sortedSamples(min([ ceil(S*ab(2)) S]));
        MCSE = (thetaB - thetaA)/2;
    % ------------------------------------------------------------------ %
    otherwise
        error('not coded yet.')
end

end