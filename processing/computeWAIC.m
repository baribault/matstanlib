function WAIC = computeWAIC(log_lik)
%COMPUTEWAIC computes WAIC from log likelihood samples. 
% 
% this function computes the Wantabe-Akaike Information Criterion (also
% known as the Widely-Applicable Information Criterion) using what is
% commonly called "pWAIC2" to compute the effetive number of parameters.
% 
% WAIC = CALCULATEWAIC(LOG_LIK)
%   this function accepts an [nIterations nChains <data dimensions>]-sized
%   matrix of log likelihood values (LOG_LIK).
%   a single value, WAIC, is computed and returned. 
% 
% 
% References: Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
%             Wantabe (2010). Journal of Machine Learning Research.
% 
% 
% (c) beth baribault 2019 ---                                 > matstanlib       

%% check inputs
if ~isnumeric(log_lik)
    error('log_lik input must be numeric.')
elseif ismatrix(log_lik)
    error(['log_lik must have at least three dimensions (i.e., should ' ...
        'be a [nIterations nChains <data dimensions>]-sized matrix).'])
end

%% WAIC
%convert [nIterations nChains <data dimensions>]
%     to [nIterations nChains nDataPoints]
log_lik = reshape(log_lik,size(log_lik,1),size(log_lik,2),[]);
%convert [nIterations nChains nDataPoints]
%     to [nSamples            nDataPoints]
log_lik = reshape(log_lik,[],size(log_lik,3));

%compute WAIC
lppd = log(mean(exp(log_lik),1));   %mean over samples
lppd = sum(lppd);                   %sum over data points
p_WAIC = var(log_lik);              %sample variance
p_WAIC = sum(p_WAIC);               %sum over data points
elppd_WAIC = lppd - p_WAIC;         %expected log pointwise predictive density
WAIC = -2*elppd_WAIC;               %convert to deviance scale

end