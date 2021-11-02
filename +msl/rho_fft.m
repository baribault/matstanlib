function [rho_hat,lags] = rho_fft(chains)
%RHO_FFT computes the combined autocorrelation across chains via FFT.
% 
% [RHO_HAT, LAGS] = RHO_FFT(CHAINS)
%   computes the combined autocorrelation across chains, RHO_HAT, 
%   for non-negative LAGS (i.e., 0:nIterations-1) from CHAINS, 
%   a [nIterations nChains]-sized matrix of posterior samples. 
%   
%   to compute the autocorrelation within each chain, the fast fourier
%   transform is used (for greater efficiency and computational accuracy). 
% 
% (c) beth baribault 2021 ---                                 > matstanlib

import msl.*

%%
N = size(chains,1); %number of (split) iterations
M = size(chains,2); %number of (split) chains

chainMeans = mean(chains);
overallMean = mean(chainMeans);
chainVariances = var(chains,0); %with N-1
betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
withinChainsVar = 1/M * sum(chainVariances);                      %W
overestPostVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+

%autocorrelation by lag, across chains
varXacf = zeros([N 1]);
for m = 1:M
    [acorr,lags] = acf(chains(:,m));
    acorr = acorr/acorr(1); %normalize ???
    varXacf = varXacf + chainVariances(m)*acorr;
end
rho_hat = 1 - (withinChainsVar - (1/M)*varXacf)/overestPostVar;

end