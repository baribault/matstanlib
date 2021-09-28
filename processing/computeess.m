function ESS = computeess(chains,method,alpha)
%COMPUTEESS computes the effective sample size (ESS) from posterior samples.  
% 
% ESS = COMPUTEESS(CHAINS,TYPE)
% ESS = COMPUTEESS(CHAINS,'quantile',ALPHA)
% ESS = COMPUTEESS(CHAINS,'local',ALPHA)
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
%     * 'local'     >>  ESS for a particular quantile range, ALPHA
% 
%   if an above TYPE is marked with *, then an additional input, ALPHA, is
%   also required.  
%   if TYPE is 'quantile', then ALPHA is a **single number**, the quantile
%   for which to compute ESS, given as a proportion.  
%   if TYPE is 'local', then ALPHA must have **two elements**, the minimum
%   and maximum quantiles bounding the local range for which to compute
%   ESS, both given as a proportion.  the first element must be lower than
%   the second. 
%   in either case, ALPHA must be (a) numeric value(s) between 0 and 1.  
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
validMethods = {'bulk','tail','mean','median','sd','quantile','local', ...
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
switch method
    case 'quantile'
        nAlpha = 1;
        noAlphaErrMsg = ...
            ['if ESS type is ''quantile'', then alpha (the ' ...
            'requested ESS quantile, as a proportion), is required.'];
    case 'local'
        nAlpha = 2;
        noAlphaErrMsg = ...
            ['if ESS type is ''local'', then alpha (a two-element vector ' ...
            'representing the inclusive minimum and maximum quantiles) ' ...
            ' is required.'];
end
if ismember(method,{'quantile','local'})
    if nargin < 3
        error(noAlphaErrMsg)
    elseif ~isnumeric(alpha) || ~isequal(numel(alpha),nAlpha)
        if nAlpha==1, error('alpha must be a single numeric value.')
        else,         error('alpha must be a two-element numeric vector.')
        end 
    elseif any(alpha <= 0 | alpha >= 1)
        error('alpha values must be in the interval (0,1).')
    elseif nAlpha==2 && alpha(1)>=alpha(2)
        error(['if alpha has two elements, the first element ' ...
            'cannot be greater than or equal to the second.'])
    end
end

%% compute ESS

switch method
    % ------------------------------------------------------------------ %
    case 'BDA2'
        ESS = msl.ess_BDA2(chains);
    % ------------------------------------------------------------------ %
    case 'BDA3'
        ESS = msl.ess_BDA3(msl.splitchains(chains));
        fprintf('BDA3      = %.5f\n',ESS)
        
        ESS = msl.ess_core(msl.splitchains(chains));
        fprintf('CORE      = %.5f\n',ESS)
        
        ESS = msl.ess_tmp(msl.splitchains(chains));
        fprintf('TMP       = %.5f [TEST]\n',ESS)
        
    % ------------------------------------------------------------------ %
    case 'bulk'
        ESS = msl.ess_core(msl.ranknorm(msl.splitchains(chains)));
    case 'tail'
        ESS_q05 = msl.ess_core(msl.ranknorm(msl.splitchains(chains)),0.95);
        ESS_q95 = msl.ess_core(msl.ranknorm(msl.splitchains(chains)),0.05);
        ESS = min([ESS_q05 ESS_q95]);
    % ------------------------------------------------------------------ %
    case 'mean'
        ESS = msl.ess_core(msl.splitchains(chains));
    case 'sd'
        ESS_m = msl.ess_core(msl.splitchains(chains));
        ESS_s = msl.ess_core(msl.splitchains(chains.^2));
        ESS = min([ESS_m ESS_s]);
    % ------------------------------------------------------------------ %
    case 'median'
        ESS = msl.ess_core(msl.splitchains(chains),0.50);
    case 'mad'
        chains = abs(chains - median(chains(:)));
        Ichains = chains <= median(chains(:));
        ESS = msl.ess_core(msl.ranknorm(msl.splitchains(Ichains)));
    case 'quantile'
        ESS = msl.ess_core(msl.splitchains(chains),alpha);
    case 'local'
        ESS = msl.ess_core(msl.splitchains(chains),alpha);
    % ------------------------------------------------------------------ %
    otherwise
        error('not coded yet.')
end

%enforce maximum for ESS estimate
if ~ismember(method,{'BDA2','BDA3'})
    S = numel(chains);  %total number of samples
    ESS = min(ESS,S*log10(S));
end

end