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
S = numel(chains);  %total number of samples

%combined autocorrelation (across chains)
[rho_hat,lags] = rho_fft(chains);

%compute P, the combined autocorrelation over pairs
T = floor(length(rho_hat)/2);
P_hat = NaN([T 1]);
for t = 0:(T-1)
    t1 = 2*t;
    t2 = 2*t + 1;
    P_hat(t+1) = rho_hat(t1+1) + rho_hat(t2+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;hold on;
% plot(lags,zeros(size(lags)),'k:')
% plot(lags,rho_hat)
% plot(0:length(P_hat)-1,P_hat)
% c = find(P_hat<=0,1);
% plot(c-1,P_hat(c),'ok','markersize',12)
% xlim([0 20])
% legend('0','rho hat','P hat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%truncate P such that all values of P are positive
idx1stNegAutocorrPair = find(P_hat<=0,1);
if ~isempty(idx1stNegAutocorrPair)
    P_hat = P_hat(1:idx1stNegAutocorrPair-1);
end

%compute tau_hat, ESS
tau_hat = -1 + 2*sum(P_hat);
ESS = S/tau_hat;

% %enforce minimum for ESS estimate ???
% ESS = max(ESS,M);

end