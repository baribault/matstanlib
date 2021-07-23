function rtable = rhattable(samples,parameterRequest)
%RHATTABLE calculates convergence metrics based on MCMC samples.
% 
% RTABLE = RHATTABLE(SAMPLES)
%   for all instances of all parameters found in SAMPLES, metrics 
%   useful for establishing convergence based on the mcmc samples alone 
%   are calculated.
%   a table, RTABLE, with the following variables is returned:
%     - parameter:    parameter instance names
%     - meanval:      the mean of the posterior samples
%     - rhat:         split rhat (potential scale reduction factor)
%     - neff:         effective sample size
%     - neff_ratio:   effective sample size ratio [0--1]
% 
% RTABLE = RHATTABLE(SAMPLES,PARAMETERREQUEST)
%   generates a RTABLE output for the requested parameters only.
%   PARAMETERREQUEST must be a valid parameter name string or cell of 
%   valid parameter name strings. instance names may be included.
% 
% NOTE: at time of writing (9/2019), this function's output matches the 
% summary output from Stan, as both use the calculations specified in 
% BDA3.
% 
% Reference:  Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also EXTRACTSAMPLES, INTERPRETDIAGNOSTICS
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% check inputs
if nargin > 2
    error('too many inputs.')
end

%parse inputs to determine parameter(s)
if nargin == 1
    parameterRequest = fieldnames(samples);
elseif nargin == 2
    if ischar(parameterRequest)
        parameterRequest = {parameterRequest};
    elseif ~iscell(parameterRequest) || ~all(cellfun(@ischar,parameterRequest))
        error(['parameter names must be input as a string or ' ...
               'cell of strings.'])
    end
end

%% calculate gelman-rubin statistics
%create a list of parameter instances
parameters = getparaminstances(parameterRequest, ...
    fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));

meanval = NaN(size(parameters));
rhat = NaN(size(parameters));
neff = NaN(size(parameters));
neff_ratio = NaN(size(parameters));
isInstance = cellfun(@(x) any(x=='['),parameters);
for p = 1:length(parameters)
    %account for parameters vs. parameter instances
    if isInstance(p)
        [parameter,ind] = str2ind(parameters{p});
        chains = samples.(parameter)(:,:,ind{:});
    else
        parameter = parameters{p};
        chains = samples.(parameter);
    end
    
    %mean value
    meanval(p) = mean(chains(:));
    
    %split chains
    halfN = floor(size(chains,1)/2);
    splitChains = [chains(1:halfN,:) chains(halfN+1:halfN*2,:)];
    
    %calculate split rhat
    N = size(splitChains,1); %number of iterations
    M = size(splitChains,2); %number of chains
    chainMeans = mean(splitChains,1);
    overallMean = mean(splitChains(:));
    chainVariances = var(splitChains,1);
    betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);
    withinChainsVar = 1/M * sum(chainVariances);
    pooledVar = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar; %vhat
    rhat(p) = sqrt(pooledVar/withinChainsVar);
    
    %calculate the number of effective samples
    rho_hat = acf_hat(splitChains,pooledVar);
    neff_hat = M*N/(1+2*sum(rho_hat));
    neff(p) = floor(min(M*N,neff_hat));
    neff_ratio(p) = neff(p)/(M*N);
end

%collect output in a table
parameter = parameters; %quick rename!
rtable = table(parameter,meanval,rhat,neff,neff_ratio);

end

%------------------------------------------------------------------------%
function rho_hat = acf_hat(splitChains,pooledVar)
%ACF_HAT calculates the autocorrelation via variogram for ESS calculation.
    %extract dimensions
    N = size(splitChains,1);
    M = size(splitChains,2);
    maxLag = N - 1;
    %calculate variogram & autocorrelation
    variogram = NaN([maxLag 1]);
    rho_hat = NaN([maxLag 1]);
    keepCalculating = true; endSum = false; lag = 0; 
    while keepCalculating
        lag = lag + 1;
        sumArg = (splitChains(1+lag:end,:) - splitChains(1:end-lag,:)).^2;
        variogram(lag) = 1/(M*(N-lag))*sum(sum(sumArg));
        rho_hat(lag) = 1 - variogram(lag)/(2*pooledVar);
        %criterion for ending the partial sum
        if lag > 2 && mod(lag,2) == 1
            endSum = (rho_hat(lag) + rho_hat(lag-1)) < 0;
        end
        %criteria for ending autocorrelation calculation
        if endSum || lag == maxLag
            keepCalculating = false;
        end
    end
    rho_hat(isnan(rho_hat)) = [];
end