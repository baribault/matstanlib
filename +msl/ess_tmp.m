function ESS = ess_core(chains,proportion)
%ESS_CORE applies the base effective sample size formula to posterior samples.
% 
% ESS = ESS_CORE(CHAINS)
%   computes the effective sample size, ESS, for CHAINS, a [nIterations
%   nChains]-sized matrix of posterior samples, based on the estimated
%   combined autocorrelation function.  
%   
%   if CHAINS has folded & split chains,  then ESS is bulk ESS.
%   if CHAINS has split chains of indicators relative to a quantile, 
%                                         then ESS is quantile ESS.
%   if CHAINS has split chains,           then ESS is BDA3 n_eff.
%   if 'ess_BDA2.m' is used instead,      then ESS is BDA2 n_eff.
% 
% References: Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% 
% See also COMPUTEESS, ESS_BDA2, MCMCTABLE
% 
% (c) beth baribault 2021 ---                                 > matstanlib 

import msl.*

%% check inputs
if nargin < 1, error('a chains matrix is required.'); end
%chains
if ~isnumeric(chains) || ~ismatrix(chains)
    error('sole input must be chains, a matrix of posterior samples.')
end

%proportion
if nargin < 2 || isempty(proportion)
    proportion = NaN; %default
elseif isscalar(proportion) && isnumeric(proportion)
    if proportion <= 0 || proportion >= 1
        error('proportion must be in the interval (0,1).')
    end
else
    error('proportion must be a single numeric value.')
end

%% if ESS is relative to a quantile, convert to indicators
if ~isnan(proportion)
    chains = chains <= quantile(chains(:),proportion);
    chains = double(chains); %not f*cking logical type, christ
end

%% compute the effective sample size

N = size(chains,1); %number of (split) iterations
M = size(chains,2); %number of (split) chains
S = numel(chains);  %total number of samples

chainMeans = mean(chains);
overallMean = mean(chainMeans);
chainVariances = var(chains,0); %[with 1/(N-1)]
betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);    %B
withinChainsVar = 1/M * sum(chainVariances);                      %W
overestPostVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;  %var+

% %autocorrelation by lag, across chains
% varXacf = zeros([N 1]);
% for m = 1:M
%     [acorr,lags] = acf(chains(:,m));
%     acorr = acorr/acorr(1); %normalize ???
%     varXacf = varXacf + chainVariances(m)*acorr;
% end
% rho_hat = 1 - (withinChainsVar - (1/M)*varXacf)/overestPostVar;
%autocorrelation by lag, across chains
%% 
acov = NaN([N 4]);
for m = 1:M
    [acorr,~] = acf_(chains(:,m));
    acov(:,m) = acorr;
    acov(:,m) = acorr/acorr(1) * var(chains(:,m),0) * (N-1)/N;
end
% acov(1,:)

mean_var = mean(acov(1,:))*N/(N - 1);
% mean_var - mean(chainVariances) %basically the same

var_plus = mean_var*(N - 1)/N + var(chainMeans);
% var_plus - overestPostVar %basically the same

% mean_var = mean(chainVariances);
% var_plus = overestPostVar;

%% Geyer's initial positive sequence
rho_hat_t = zeros([N 1]);
t = 0;
% rho_hat_even = 1 - (mean_var - mean(acov(t+1,:)))/var_plus
rho_hat_even = 1;
rho_hat_t(t+1) = rho_hat_even;

rho_hat_odd  = 1 - (mean_var - mean(acov(t+2,:)))/var_plus;
rho_hat_t(t+2) = rho_hat_odd;

t = 0;
while t < size(acov,1) - 5 && ~isnan(rho_hat_even + rho_hat_odd) && ...
        rho_hat_even + rho_hat_odd > 0
    t = t + 2;
    rho_hat_even = 1 - (mean_var - mean(acov(t+1,:))) / var_plus;
    rho_hat_odd  = 1 - (mean_var - mean(acov(t+2,:))) / var_plus;
    if rho_hat_even + rho_hat_odd >= 0
        rho_hat_t(t+1) = rho_hat_even;
        rho_hat_t(t+2) = rho_hat_odd;
    end
end
max_t = t;

% this is used in the improved estimate
if rho_hat_even > 0
    rho_hat_t(max_t + 1) = rho_hat_even;
end

%% Geyer's initial monotone sequence
t = 0;
while (t <= max_t - 4)
    t = t + 2;
    if rho_hat_t(t+1) + rho_hat_t(t+2) > rho_hat_t(t-1) + rho_hat_t(t)
        rho_hat_t(t+1) = (rho_hat_t(t-1) + rho_hat_t(t)) / 2;
        rho_hat_t(t+2) = rho_hat_t(t+1);
    end
end

%% Geyer's truncated estimate
% tau_hat = -1 + 2*sum(rho_hat_t(1:max_t))
%Improved estimate reduces variance in antithetic case
tau_hat = -1 + 2*sum(rho_hat_t(1:max_t)) + rho_hat_t(max_t+1);

%Safety check for negative values and with max ess equal to ess*log10(ess)
tau_hat = max(tau_hat, 1/log10(S));

ESS = S / tau_hat;

end