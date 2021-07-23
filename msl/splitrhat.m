function PSRF = splitrhat(chains,alreadySplit)
%SPLITRHAT calculates split rhat based on a matrix of MCMC samples.  
%   
% PSRF = SPLITRHAT(CHAINS)
%   calculates the PSRF, the potential scale reduction factor (i.e., Rhat),
%   for CHAINS, a [nIterations nChains]-sized matrix of posterior samples.
%   the split Rhat value, PSRF, is returned. 
% 
% PSRF = SPLITRHAT(CHAINS,ALREADYSPLIT)
%   this syntax may be used to prevent CHAINS that have already been 
%   split being split a second time (or to force a non-split Rhat
%   calculation).
%   if ALREADYSPLIT is false or 0 , then CHAINS will be split.
%   if ALREADYSPLIT is true or 1, then CHAINS will not be split.
%   the default value of ALREADYSPLIT is false.
% 
% References: Gelman & Rubin (1992). Statistical Science.
%             Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also MCMCTABLE
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% check inputs
if nargin < 1, error('a chains matrix is required.'); end
%chains
if ~isnumeric(chains) || ~ismatrix(chains)
    error('')
end

%alreadySplit
if nargin == 1
    alreadySplit = false;
elseif nargin == 2
    if ~islogical(alreadySplit) && ~ismember(alreadySplit,[0 1])
        error('the alreadySplit input must be a logical value or 0 or 1.')
    end
end

%% compute PSRF
%split chains
if ~alreadySplit
    splitChains = splitchains(chains);
end

%calculate potential scale reduction factor from split chains
N = size(splitChains,1); %number of (split) iterations
M = size(splitChains,2); %number of (split) chains
chainMeans = mean(splitChains);
overallMean = mean(splitChains(:));
chainVariances = var(splitChains,0);
betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
withinChainsVar = 1/M * sum(chainVariances);                      %W
margPostVarEst = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+
PSRF = sqrt(margPostVarEst/withinChainsVar);

end