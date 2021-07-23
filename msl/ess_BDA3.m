function ESS = ess_BDA3(chains)
%ESS_BDA3 computes the effective sample size as per BDA3. 
% 
% *** THIS COMPUTATION SHOULD BE USED FOR BACKWARDS COMPATABILITY ONLY! ***
% 
% N_EFF = ESS_BDA3(CHAINS)
%   the matrix of posterior samples, CHAINS, is used to compute the
%   effective sample size, N_EFF, as per BDA3.  
% 
% Reference:  Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% (c) beth baribault 2021 ---                                 > matstanlib 

%% check inputs
if nargin < 1, error('a chains matrix is required.'); end
%chains
if ~isnumeric(chains) || ~ismatrix(chains)
    error('sole input must be chains, a matrix of posterior samples.')
end

%% compute the effective sample size
N = size(chains,1); %number of (split) iterations
S = numel(chains);  %total number of samples

%compute the "variogram"
[rho_hat,lags] = rho_variogram(chains);

%compute P, the combined autocorrelation over pairs
startIdx = 3:2:N-1;
n = 0;
notNegative = true;
while notNegative
    n = n + 1;
    ind1 = startIdx(n); ind2 = startIdx(n)+1;
    if rho_hat(ind1) + rho_hat(ind2) <=0
        T = lags(startIdx(n)-1);
        notNegative = false;
    end
end

%truncate rho_hat to "first odd integer" ...
rho_hat = rho_hat(ismember(lags,1:T));

%compute tau_hat, ESS
tau_hat = 1 + 2*sum(rho_hat);
ESS = S/tau_hat;

end

%% --------------------------------------------------------------------- %%
function [rho_hat,lags] = rho_variogram(chains)
%RHO_VARIOGRAM computes the combined autocorrelation across chains via variogram.
% 
% *** THIS COMPUTATION SHOULD BE USED FOR BACKWARDS COMPATABILITY ONLY! ***
% 
% [RHO_HAT, LAGS] = RHO_VARIOGRAM(CHAINS)
%   computes the combined autocorrelation across chains, RHO_HAT, 
%   for non-negative LAGS (i.e., 0:nIterations-1) from CHAINS, 
%   a [nIterations nChains]-sized matrix of posterior samples.  
%   
%   this function uses the variogram method of BDA3. 
% 
% (c) beth baribault 2021 ---                                 > matstanlib
    
    N = size(chains,1); %number of (split) iterations
    M = size(chains,2); %number of (split) chains
    
    %compute var+
    chainMeans = mean(chains);
    overallMean = mean(chains(:));
    chainVariances = var(chains,0); %normalize by N-1
    betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
    withinChainsVar = 1/M * sum(chainVariances);                      %W
    overestPostVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+
    
    %compute "variogram"
    lags = (0:(N-1))';
    vario = NaN([N 1]);
    for n = 1:N
        t = lags(n);
        doublesum = 0;
        for j = 1:M
            for i = (t+1):N
                doublesum = doublesum + (chains(i,j) - chains(i-t,j)).^2;
            end
        end
        vario(n) = 1/(M*(N-t))*doublesum;
    end
    %... and "invert the formula" ...
    rho_hat = 1 - vario/(2*overestPostVar);
end