function rhat = rhat(chains,varargin)
%RHAT computes the Rhat diagnostic from posterior samples.  
%   
% RHAT = RHAT(CHAINS)
%   computes the Rhat diagnostic from a [nIterations nChains]-sized matrix
%   of posterior samples, CHAINS.  
%   the Rhat value, RHAT, is returned.
% 
% 
% there are two optional inputs, which may be given in any order.  
% optional inputs are parsed by type.  
% 
% RHAT = RHAT(...,METHOD)
%   by default, the current method for computing Rhat (as per Vehtari et
%   al., 2020) is used.  for backwards compatability, we include the
%   previous standard method for computing Rhat (as per BDA3).  
%   as such, valid options for METHOD include:
%       'vehtari'   >>  the current method
%       'BDA3'      >>  the method described in BDA3
% 
% RHAT = RHAT(...,ALREADYSPLIT)
%   this optional input is useful to accomodate a CHAINS input that has
%   already been split (or to force a non-split computation). 
%   if ALREADYSPLIT is false or 0 , then CHAINS will be split.
%   if ALREADYSPLIT is true or 1, then CHAINS will not be split.
%   the default value of ALREADYSPLIT is false.
% 
% 
% Reference for current Rhat formulae:
%             Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% Reference for 'BDA3' formulae:
%             Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                 Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also MCMCTABLE
% 
% (c) beth baribault 2021 ---                                 > matstanlib

import msl.*

%% check inputs
if ~isnumeric(chains) && ~ismatrix(chains)
    error(['the first input must be chains, a [nIterations nChains]-sized ' ...
        'matrix of posterior samples'])
end

%%% optional inputs %%%                
method = 'vehtari';     foundMethod = false;        validMethods = {'vehtari','BDA3'};
alreadySplit = false;   foundAlreadySplit = false;      

for v = 1:length(varargin)
    %mode
    if ischar(varargin{v})
        if ismember(varargin{v},validMethods)
            if ~foundMethod,  method = varargin{v};
            else,           error('only one mode input is permitted.')
            end
        else, error(['computation method ''%s'' is not valid.  valid ' ...
                'inputs include: ''%s''.'],strjoin(validMethods,''', '''))
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

%% compute rhat

% N = size(chains,1); %number of (split) iterations
% M = size(chains,2); %number of (split) chains
% S = numel(chains);  %total number of samples

switch method
    case 'BDA2'
        rhat = splitrhat(chains);
    case 'BDA3'
        chains = splitchains(chains);
        rhat = splitrhat(chains);
    case 'vehtari'
        chains = splitchains(chains);
        rhat_bulk = splitrhat(ranknorm(chains));
        chains = chains - median(chains(:));   %folded draws
        rhat_tail = splitrhat(ranknorm(chains));
        rhat = max(rhat_bulk,rhat_tail);
end

end

% ---------------------------------------------------------------------- %
function chains = splitchains(chains)
    halfN = floor(size(chains,1)/2);
    chains = [chains(1:halfN,:) chains(end-halfN+1:end,:)];
end

% ---------------------------------------------------------------------- %
function psrf = splitrhat(chains)
    %as per BDA3 and Vehtari et al. (2020)
    N = size(chains,1); %number of (split) iterations
    M = size(chains,2); %number of (split) chains
    chainMeans = mean(chains);
    overallMean = mean(chainMeans);
    chainVariances = var(chains,0); %[with 1/(N-1)]
    betweenChainsVar = N/(M-1) * sum((chainMeans-overallMean).^2);
    withinChainsVar = 1/M * sum(chainVariances);
    margPostVarEst = (N-1)/N*withinChainsVar + 1/N*betweenChainsVar;
    psrf = sqrt(margPostVarEst/withinChainsVar);
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