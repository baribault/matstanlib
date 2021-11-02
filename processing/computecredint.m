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
%   *** NOTE: if the requested % CI is so low that there are not enough
%   *** samples to offer an interval (i.e., PROPORTION is < 2/nSamples),
%   *** then CI will be [NaN NaN]. 
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

% %% check there are enough samples to actually compute an *interval*
% nSamples = numel(chains);
% if proportion < 2/nSamples
%     error(['there are not enough samples to compute a %g%% credible interval!' ...
%         '(%i samples are available, but at least %i would be needed.)'], ...
%         proportion*100,nSamples,ceil(2/nSamples*100))
% end

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
        nInCI = max(floor(nSamples*proportion),1);
        nIntervals = nSamples - nInCI + 1;
        %which of all possible intervals is the shortest?
        orderedSamples = sort(chains,'ascend');
        width = NaN([nIntervals 1]);
        for n = 1:nIntervals 
            width(n) = diff([orderedSamples(n) orderedSamples(n + nInCI - 1)]);
        end
        minWidth = min(width);
        %cope with one vs. many possible shortest intervals
        startIndices = find(width==minWidth);
        if isscalar(startIndices)
            %if there is a single shortest HDI
            start = startIndices;
        else
            %if multiple shortest HDI, pick the one in the middle
            [~,start] = min(abs(startIndices+minWidth/2 - ceil(nSamples/2)));
        end
        %return HDI
        CI = [orderedSamples(start) orderedSamples(start + nInCI -1)];
end

end