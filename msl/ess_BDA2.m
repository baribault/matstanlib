function n_eff = ess_BDA2(chains)
%ESS_BDA2 computes the effective sample size as per BDA2. 
% 
% *** THIS COMPUTATION SHOULD BE USED FOR BACKWARDS COMPATABILITY ONLY! ***
% 
% N_EFF = ESS_BDA2(CHAINS)
%   the matrix of posterior samples, CHAINS, is used to compute the
%   effective sample size, N_EFF, as per BDA2.  
% 
% Reference:  Gelman, Carlin, Stern, & Rubin (2003). Bayesian Data
%                 Analysis, 2nd edition. Chapman & Hall.
% 
% (c) beth baribault 2021 ---                                 > matstanlib

import msl.*

%%
N = size(chains,1); %number of iterations
M = size(chains,2); %number of chains

chainMeans = mean(chains);
overallMean = mean(chains(:));
chainVariances = var(chains,0); %normalize by N-1
betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
withinChainsVar = 1/M * sum(chainVariances);                      %W
overestPostVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+

n_eff = M*N*(overestPostVar/betweenChainsVar);
n_eff = min(M*N,n_eff); %maximum reportable value is number of samples

end