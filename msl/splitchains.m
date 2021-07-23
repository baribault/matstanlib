function chains = splitchains(chains)
% SPLITCHAINS splits a matrix of posterior samples by chains. 
% 
% CHAINS = SPLITCHAINS(CHAINS)
%   CHAINS, the [nIterations nChains]-sized matrix of posterior samples, is
%   split such that each chain (i.e., column) is divided in two, resulting
%   in a [floor(nIterations/2) 2*nChains]-sized matrix. 
%   
%   NOTE: if there are an odd number of iterations, then the last sample in
%   each chain is discarded (as in the R packages cmdstanr & posterior).  
%   a different approach is used in the Python library Arviz, which
%   discards the middle sample in each chain. 
% 
% (c) beth baribault 2021 ---                                 > matstanlib

N = floor(size(chains,1)/2); %number of (split) iterations

% chains = [chains(end-2*N+1:end-N,:) chains(end-N+1:end,:)]; %skip first sample
chains = [chains(1:N,:) chains(N+1:2*N,:)];       %skip last sample, as in R (posterior)
% chains = [chains(1:N,:) chains(end-N+1:end,:)]; %skip middle sample, as in Python (Arviz)

end