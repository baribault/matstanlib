function ESS = computeess(chains,method,alpha)
%COMPUTEESS computes the effective sample size (ESS) from posterior samples.  
%   
% ESS = COMPUTEESS(CHAINS,TYPE)
%   this function computes ESS, the estimated effective sample size, from 
%   CHAINS, a [nIterations nChains]-sized matrix of posterior samples.
%   
%   TYPE determines which type of ESS estimate should be computed.
%   many current estimate types are available:
%       'bulk'      >>  bulk ESS (for the more central parameter values)
%       'tail'      >>  tail ESS (for the more extreme parameter values)
%       'mean'      >>  ESS for the posterior mean
%       'median'    >>  ESS for the posterior median
%       'sd'        >>  ESS for the posterior standard deviation
%     * 'quantile'  >>  ESS for a particular quantile, ALPHA
%   
% ESS = COMPUTEESS(CHAINS,'quantile',ALPHA)
%   if TYPE is 'quantile', then an additional input, ALPHA, is required.
%   ALPHA is the quantile for which to compute ESS, given as a proportion. 
%   as such, ALPHA must be a scalar numeric value between 0 and 1.  
%   
%   for all of the above, current best computational practices are used. 
% 
% 
% ESS = COMPUTEESS(CHAINS,METHOD)
%   for backwards compatability, previously used methods of computing ESS 
%   are also available.  
%   to use an older computation METHOD, give one of the following strings
%   as an alterenative input to TYPE:
%       'BDA3'      >>  N_eff, as per BDA3
%       'BDA2'      >>  N_eff, as per BDA2
% 
% 
% Reference for current ESS computation methods:
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% Reference for 'BDA3' computation method:
%             Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% Reference for 'BDA2' computation method:
%             Gelman, Carlin, Stern, & Rubin (2003). Bayesian Data
%                 Analysis, 2nd edition. Chapman & Hall.
% 
% 
% See also MCMCTABLE, COMPUTERHAT, COMPUTEMCSE
% 
% (c) beth baribault 2021 ---                                 > matstanlib

%% check inputs
%chains
if nargin < 1, error('a chains matrix is required.'); end
if ~isnumeric(chains) && ~ismatrix(chains)
    error(['the first input must be chains, a [nIterations nChains]-sized ' ...
        'matrix of posterior samples'])
end

%method
validMethods = {'bulk','tail','mean','median','sd','quantile', ...
                'BDA3','BDA2'};
if nargin < 2
    %type is REQUIRED!
    error(['an ESS type or computation method is required. ' ...
        'valid computation methods include: ''%s'''], ...
        strjoin(validMethods,''', '''))
elseif ischar(method)
    if ~ismember(method,validMethods)
        error(['''%s'' is not a valid ESS computation method.  \n' ...
            'valid computation methods include: ''%s'''], ...
            strjoin(validMethods,''', '''))
    end
else
    error('ESS computation method must be a string (@ischar==true).')
end

%alpha
if isequal(method,'quantile')
    if nargin < 3
        error(['if ESS type is ''quantile'', then alpha (the ' ...
            'requested ESS quantile, as a proportion), is required.']) 
    elseif ~isnumeric(alpha) || ~isscalar(alpha)
        error('alpha must be a single numeric value.')
    elseif alpha <= 0 || alpha >= 1
        error('alpha must be in the interval (0,1).')
    end
end

%% compute ESS

switch method
    % ------------------------------------------------------------------ %
    case 'BDA2'
        ESS = ess_BDA2(chains);
    % ------------------------------------------------------------------ %
    case 'BDA3'
        ESS = ess_BDA3(splitchains(chains));
        fprintf('BDA3            = %.5f\n',ESS)
        ESS = ess_core(splitchains(chains));
        fprintf('bulk split      = %.5f\n',ESS)
        ESS = ess_tmp(splitchains(chains));
        fprintf('bulk split       = %.5f [TEST]\n',ESS)
        ESS = ess_core(ranknorm(splitchains(chains)));
        fprintf('bulk split+norm = %.5f\n',ESS)
        ESS = ess_tmp(ranknorm(splitchains(chains)));
        fprintf('bulk split+norm = %.5f [TEST]\n',ESS)
        
        
        ESS = ess_tmp(splitchains(chains));
    % ------------------------------------------------------------------ %
    case 'bulk'
        ESS = ess_core(ranknorm(splitchains(chains)));
    case 'tail'
        ESS_q5 = ess_core(ranknorm(splitchains(chains)),0.95);
        ESS_q95 = ess_core(ranknorm(splitchains(chains)),0.05);
        ESS = min([ESS_q5 ESS_q95]);
    case 'median'
        ESS = ess_core(splitchains(chains),0.50);
    case 'quantile'
        ESS = ess_core(splitchains(chains),alpha);
    otherwise
        error('not coded yet.')
end

%enforce maximum for ESS estimate
if ~ismember(method,{'BDA2','BDA3'})
    S = numel(chains);  %total number of samples
    ESS = min(ESS,S*log10(S));
end

end