function ESS = ess(chains,varargin)
%ESS computes the effective sample size (ESS) from posterior samples.  
%   
% [BULK_ESS TAIL_ESS] = ESS(CHAINS)
%   computes the bulk effective sample size (BULK_ESS) and tail effective
%   sample size (TAIL_ESS) from a [nIterations nChains]-sized matrix of
%   posterior samples, CHAINS.  
%   the Rhat value, RHAT, is returned.
%   
%   by default, the current method for computing Rhat, as per 
% 
% ESS = ESS(CHAINS,ALREADYSPLIT)
%   this syntax may be used to prevent CHAINS that have already been 
%   split being split a second time (or to force a non-split computation).
%   if ALREADYSPLIT is false or 0 , then CHAINS will be split.
%   if ALREADYSPLIT is true or 1, then CHAINS will not be split.
%   the default value of ALREADYSPLIT is false.
% 
% 
% Reference for current ESS formulae:
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% Reference for 'BDA3' formulae:
%             Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also RHAT, MCSE
% 
% (c) beth baribault 2021 ---                                 > matstanlib

import msl.*

%% check inputs
if ~isnumeric(chains) && ~ismatrix(chains)
    error(['the first input must be chains, a [nIterations nChains]-sized ' ...
        'matrix of posterior samples'])
end

%%% optional inputs %%%
mode = '';              foundMode = false;          validModes = {'bulk','tail'};
method = 'vehtari';     foundMethod = false;        validMethods = {'vehtari','BDA3','BDA2'};
alreadySplit = false;   foundAlreadySplit = false;      

for v = 1:length(varargin)
    %mode
    if ischar(varargin{v})
        if ismember(varargin{v},validModes)
            if ~foundMode,  mode= varargin{v};
            else,           error('only one ESS mode input is permitted.')
            end
    %method
        elseif ismember(varargin{v},validMethods)
            if ~foundMethod,  method = varargin{v};
            else, error('only one computation method input is permitted.')
            end
        else, error(['''%s'' is not a valid computation method or ESS mode.  ' ...
                'valid inputs include: ''%s'',''%s''.'], ...
                strjoin(validMethods,''', '''),strjoin(validModes,''', '''))
        end
    %alreadySplit
    elseif (isnumeric(varargin{v}) || islogical(varargin{v})) && ...
            ismember(varargin{v},[0 1])
        if ~foundAlreadySplit,      alreadySplit = varargin{v};
        else, error('only one logical-type input (alreadySplit) is permitted.')
        end
    % ...
    else, error('optional input type not recognized.')
    end
end

if isequal(method,'vehtari') && isempty(mode)
    error('an ESS mode is required for the current ESS computations.')
end

%% compute ESS
N = size(chains,1); %number of (split) iterations
M = size(chains,2); %number of (split) chains
S = numel(chains);  %total number of samples

switch method
    case 'BDA2'
        if alreadySplit
            error(['BDA2 ESS is computed from unsplit chains, ' ...
                'but these chains have already been split.'])
        end
        ESS = ess_BDA2(chains);
        
    case 'BDA3'
        if ~alreadySplit
            chains = splitchains(chains);
        end
        ESS = ess_core(chains);
        
    case 'vehtari'
        if ~alreadySplit
            chains = splitchains(chains);
        end
        switch mode
            case 'bulk'
                ESS = ess_core(ranknorm(chains));
            case 'tail'
                ESS = NaN;
        end
        %enforce maximum for ESS estimate
        ESS = min(ESS,S*log10(S));
end


end

% ---------------------------------------------------------------------- %
function ESS = ess_core(chains)
    M = size(chains,2); %number of (split) chains
    S = numel(chains);  %total number of samples
    [rho_hat,lags] = acf(chains);
    %%
    figure;hold on;
    plot(lags,zeros(size(lags)),'k:');plot(lags,rho_hat); xlim([0 20])
    %%
    %sum over pairs
    T = floor(length(rho_hat)/2);
    P = NaN([T 1]);
    for t = 0:(T-1)
        t1 = 2*t;
        t2 = 2*t + 1;
        P(t+1) = rho_hat(t1+1) + rho_hat(t2+1);
    end
    %%
    figure;hold on;
    plot(0:length(P)-1,zeros(size(P)),'k:');plot(0:length(P)-1,P); xlim([0 20])
    %truncate such that all values of P are positive
    firstNegAutocorrPair = find(P<=0,1);
    if ~isempty(firstNegAutocorrPair)
        P = P(1:firstNegAutocorrPair-1);
    end
    %compute ESS
    tau_hat = -1 + 2*sum(P)
    ESS = S/tau_hat;
    %enforce minimum for ESS estimate
%     ESS = max(ESS,M);
end

% ---------------------------------------------------------------------- %
function [rho_hat,lags] = acf(chains)
    %%
    N = size(chains,1); %number of (split) iterations
    M = size(chains,2); %number of (split) chains
    S = numel(chains);  %total number of samples
    chainMeans = mean(chains);
    overallMean = mean(chainMeans);
    chainVariances = var(chains,0); %[with 1/(N-1)]
    betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);  %B
    withinChainsVar = 1/M * sum(chainVariances);                    %W
    margPostVarEst = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;%var+
    %autocorrelation by lag, across chains
    varXacf = zeros([2*N-1 1]);
    for m = 1:M
%         x = chains(:,m);
%         S = numel(x);       %total number of samples
%         x = x(:);           %convert from matrix to vector
%         %%% trying not to require Signal Processing Toolbox for xcorr %%%
%         %%% so cheers to stack overflow for the following few lines!  %%%
%         %# autocorrelation
%         nfft = 2^nextpow2(2*S-1);
%         acorr = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );
%         %# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
%         lags = (-(S-1):(S-1))';
%         acorr = [acorr(end-S+2:end) ; acorr(1:S)];
        [acorr,lags] = acf_fft(chains(:,m));
        varXacf = varXacf + chainVariances(m)*acorr;
    end
    rho_hat = 1 - (withinChainsVar - (1/M)*varXacf)/margPostVarEst;
    %return only positive lags
    rho_hat = rho_hat(lags>=0);
    lags    = lags(lags>=0);
end

% ---------------------------------------------------------------------- %
function [acorr,lags] = acf_fft(x)
    S = numel(x);       %total number of samples
    x = x(:);           %convert from matrix to vector
    %%% trying not to require Signal Processing Toolbox for xcorr %%%
    %%% so cheers to stack overflow for the following few lines!  %%%
    %# autocorrelation
    nfft = 2^nextpow2(2*S-1);
    acorr = ifft( fft(x,nfft) .* conj(fft(x,nfft)) );
    %# rearrange and keep values corresponding to lags: -(len-1):+(len-1)
    lags = (-(S-1):(S-1))';
    acorr = [acorr(end-S+2:end) ; acorr(1:S)];
    % %compare with MATLAB's XCORR output
    % all( (xcorr(x)-acor) < 1e-10 )
end

% ---------------------------------------------------------------------- %
function chains = splitchains(chains)
    halfN = floor(size(chains,1)/2);
    chains = [chains(1:halfN,:) chains(end-halfN+1:end,:)];
end

% ---------------------------------------------------------------------- %
function ess = ess_BDA2(chains)
    %as per BDA2 and Vehtari et al. (2020)
    N = size(chains,1); %number of (split) iterations
    M = size(chains,2); %number of (split) chains
    chainMeans = mean(chains);
    overallMean = mean(chainMeans);
    chainVariances = var(chains,0); %[with 1/(N-1)]
    betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);  %B
    withinChainsVar = 1/M * sum(chainVariances);                    %W
    margPostVarEst = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;%var+
    ess = M*N*margPostVarEst/betweenChainsVar;
end

% ---------------------------------------------------------------------- %
function z = ranknorm(x,offset)
    if nargin < 2
        %as per Vehtari et al. (2020), with reference to Blom (1958)
        offset = 3/8;
    end
    sz = size(x);       %extract size
    S = prod(sz);       %total number of samples
    x = x(:);           %convert from matrix to vector
    x = tiedrank(x);    %rank
    x = (x - offset)/(S - 2*offset + 1); %fractional offset
    z = norminv(x,0,1); %normalize
    z = reshape(z,sz);  %convert from vector to matrix
end