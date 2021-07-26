function RHAT = computerhat(chains,method)
%RHAT computes the Rhat convergence statistic from posterior samples.  
% 
% RHAT = COMPUTERHAT(CHAINS)
%   this function computes RHAT, the Rhat convergence diagnostic, from
%   CHAINS, a [nIterations nChains]-sized matrix of posterior samples.
%   
%   by default, the current best computational practice is used. 
% 
% 
% RHAT = COMPUTERHAT(CHAINS,METHOD)
%   for backwards compatability, previously used methods of computing Rhat
%   are also available.  
%   to use a specific computation METHOD, give one of the following strings
%   as an additional input: 
%       'vehtari' or 'folded' or 'folded-split' >>  folded-split Rhat
%       'BDA3'    or 'split'                    >>  split Rhat
%       'BDA2'    or 'unsplit'                  >>  Rhat
%   the default value for METHOD is 'current'.
% 
% 
% Reference for 'current' computation method:
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020, May). 
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
% See also MCMCTABLE, COMPUTEESS, COMPUTEMCSE
% 
% (c) beth baribault 2021 ---                                 > matstanlib

%% check inputs
%chains
if nargin < 1, error('a chains matrix is required.'); end
if ~isnumeric(chains) && ~ismatrix(chains)
    error(['the first input must be chains, a [nIterations nChains]-sized ' ...
        'matrix of posterior samples.'])
end

%%% optional inputs %%%
%method
validMethods = {'current', ...
                'vehtari','folded-split','folded', ...  %folded-split Rhat
                'BDA3',   'split', ...                  %split Rhat
                'BDA2',   'unsplit'};                   %Rhat
if nargin < 2
    method = 'current'; %default is current best practice
elseif ischar(method)
    if ~ismember(method,validMethods)
        error(['''%s'' is not a valid Rhat computation method.  \n' ...
            'valid computation methods include: ''%s'''], ...
            method,strjoin(validMethods,''', '''))
    end
else
    error('Rhat computation method must be a string (@ischar==true).')
end

%% compute Rhat
switch method
    case {'current','vehtari','folded-split','folded'}
        rhat_bulk = psrf(ranknorm(splitchains( ...
                         abs(chains - median(chains(:))) )));   %folded
        rhat_tail = psrf(ranknorm(splitchains(chains)));        %not folded
        RHAT = max(rhat_bulk,rhat_tail);
        
    case {'BDA3','split'}
        RHAT = psrf(splitchains(chains));
        
    case {'BDA2','unsplit'}
        RHAT = psrf(chains);
end

end