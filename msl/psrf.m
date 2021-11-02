function PSRF = psrf(chains)
%PSRF applies the base Rhat forumla to a matrix of posterior samples.  
%   
% PSRF = PSRF(CHAINS)
%   compute PSRF, the potential scale reduction factor (i.e., Rhat), for
%   CHAINS, a [nIterations nChains]-sized matrix of posterior samples. 
%   
%   if CHAINS has folded & split chains,  then PSRF is folded-split-Rhat.
%   if CHAINS has split chains,           then PSRF is BDA3 split-Rhat.
%   if CHAINS has the original chains,    then PSRF is BDA2 Rhat.
% 
% References: Gelman & Rubin (1992). Statistical Science.
%             Gelman, Carlin, Stern, & Rubin (2003). Bayesian Data
%                 Analysis, 2nd edition. Chapman & Hall.
%             Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also COMPUTERHAT, MCMCTABLE
% 
% (c) beth baribault 2021 ---                                 > matstanlib

import msl.*

%% check inputs
if nargin < 1, error('a chains matrix is required.'); end
%chains
if ~isnumeric(chains) || ~ismatrix(chains)
    error('sole input must be chains, a matrix of posterior samples.')
end

%% compute the "potential scale reduction factor"
N = size(chains,1); %number of (potentially split) iterations
M = size(chains,2); %number of (potentially split) chains
chainMeans = mean(chains);
overallMean = mean(chains(:));
chainVariances = var(chains,0); %normalize by N-1
betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
withinChainsVar = 1/M * sum(chainVariances);                      %W
overestPostVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+
PSRF = sqrt(overestPostVar/withinChainsVar);

end