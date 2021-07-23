function CI = computecredint(chains,varargin)
%COMPUTECREDINT computes a credible interval from MCMC samples. 
% 
% CI = COMPUTECREDINT(CHAINS)
%   this function computes a credible interval from CHAINS, a matrix of
%   posterior samples for a single parameter (i.e., a scalar parameter or
%   one instance of a nonscalar parameter). 
%   by default, the 95% central (i.e., equal-quantile-tailed) credible
%   interval is computed.  
%   CI, a two-element vector containing the lower and upper bound of the 
%   credible interval, is returned. 
% 
% 
% there following optional inputs may be given in any order as they are
% recognized by type: 
% 
% CI = COMPUTECREDINT(...,PROPORTION)
%   PROPORTION specifies a different % credible interval, as a propotion.
%   as such, PROPORTION must be between 0 and 1 (but cannot be 0).
%   the default value for PROPORTION is 0.95.
%   
% CI = COMPUTECREDINT(...,CIMETHOD)
%   CIMETHOD indicates how the credible interval is to be computed. 
%   valid CIMETHOD strings include:
%       'central'   >>  central (i.e., equal-tailed) interval
%       'HDI','HPD' >>  highest posterior density (i.e., shortest) interval
%   the default value for CIMETHOD is 'central'.
%   
%   NOTE: this function will NOT offer disjoint HDIs as a solution.  
% 
% 
% See also GETSAMPLESTATS, PLOTINTERVALS, PLOTRECOVERY
% 
% (c) beth baribault 2021 ---                                 > matstanlib


%% check inputs
%chains
if ~isnumeric(chains)
    error('chains must be a numeric type.')
elseif any(isnan(chains(:)))
    error('chains cannot contain NaN values.')
else
    if ~(ismatrix(chains) || isvector(chains))
        error(['chains must be a 2D matrix of MCMC samples for ' ...
            'a *single* parameter.'])
    end
end

%%% optional %%%
proportion = 0.95;          foundProportion = false;
validCIMethods = {'central','hdi','hpd'};
CImethod = 'central';       foundCImethod = false;

for v = 1:length(varargin)
    if isempty(varargin{v})
        %(do nothing)
    %proportion
    elseif isnumeric(varargin{v}) && isscalar(varargin{v})
        proportion = varargin{v};
        if foundProportion
        	error('only one scalar numeric input (proportion) is permitted.')
        elseif proportion <= 0 || proportion > 1
            error('proportion must be in the interval (0,1].')
        end
        foundProportion = true;
    %CImethod
    elseif ischar(varargin{v})
        CImethod = lower(varargin{v});
        if foundCImethod
        	error('only one string input (CImethod) is permitted.')
        elseif ~ismember(CImethod,validCIMethods)
            error(['CImethod input ''%'' is not valid. valid options ' ...
                'include: ''%s''.'],CImethod,strjoin(validCIMethods,''', '''))
        end
        foundCImethod = true;
    else
        error('unrecognized optional input type.')
    end
end

%% compute credible interval
chains = chains(:);

switch CImethod
    % ------------------------------------------------------------------ %
    case {'central'}
        q = [(1-proportion)/2, 1 - (1-proportion)/2];
        CI = [quantile(chains,q(1)) quantile(chains,q(2))];
    % ------------------------------------------------------------------ %
    case {'hdi','hpd'}
        %how many samples required for this % CI?
        nSamples = numel(chains);
        nInCI = floor(nSamples*proportion); 
        nIntervals = nSamples - nInCI;
        %which of all possible intervals is the shortest?
        orderedSamples = sort(chains,'ascend');
        width = NaN([nIntervals 1]);
        for n = 1:nIntervals 
            width(n) = diff([orderedSamples(n) orderedSamples(n + nInCI)]);
        end
        [~,start] = min(width);
        CI = [orderedSamples(start) orderedSamples(start + nInCI)];
end

end