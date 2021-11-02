function [BF10,BF01] = savagedickey(priorChains,posteriorChains, ...
    compareAt,varargin)    
%SAVAGEDICKEY computes a Savage-Dickey Bayes factor.
% 
% this function computes a Savage-Dickey Bayes factor, which is built from
% the assumption that two nested models (represeting hypotheses H1 vs. H0)
% are being compared. 
% 
% Reference: Wagenmakers, E. J., Lodewyckx, T., Kuriyal, H., & 
%            Grasman, R. (2010). Bayesian hypothesis testing for 
%            psychologists: A tutorial on the Savage–Dickey method. 
%            Cognitive psychology, 60(3), 158-189.
% 
% the null hypothesis (H0) is that the parameter of interest is equal to a
% given value  (e.g., delta = 0).  
% the alternative hypothesis (H1) that the parameter is different from that
% value (delta ~ F).  
% the computed Bayes factor (BF_10) quantifies the relative evidenital
% support for H1 vs H0.  
% 
% NOTE: as a side effect, a figure is generated that presents the 
% smoothed prior and posterior densities and highlights their height at 
% the comparison value. 
% 
% 
% [BF10,BF01] = SAVAGEDICKEY(PRIORCHAINS,POSTERIORCHAINS,COMPAREAT)
%   PRIORCHAINS and POSTERIORCHAINS each must be a matrix of posterior
%   samples of size [nIterations nChains].  
%   COMPAREAT is the value at which the prior and posterior densities
%   will be compared.  
%   the Bayes factor comparing the alternative to the null (BF10) and 
%   the Bayes factor comparing the null to the alternative (BF01) are 
%   returned. 
% 
% 
% optional inputs may be given as name-value pairs:
% SAVAGEDICKEY(... ,'name1',VALUE1,'name1',VALUE2,...)
% pairs may be given in any order.  property names are case insensitive.
%   
%   'xi'
%       XI offers a different domain for the kernel density smoothing.  
%       XI must be a vector of numeric values (but see ksdensity.m for
%       particular details).
%   
%   'support'
%       SUPPORT offers a different support for the kernel density
%       smoothing.  SUPPORT may be 'unbounded', 'positive', or a
%       two-element vector representing the minimum & maximum value.  
% 
% 
% NOTE: this function requires ksdensity.m (part of MATLAB's Statistics 
% Toolbox) for its calculations.
% 
% 
% See also EXTRACTSAMPLES, KSDENSITY
% 
% (c) beth baribault 2019 ---                                 > matstanlib

msl.options

%% parse inputs
%value at which we will compare the prior vs posterior for delta
if nargin == 2
    %default value for compareAt
    compareAt = 0;
elseif nargin > 2
    %check the given compareAt input
    if ~isnumeric(compareAt) || ~isscalar(compareAt)
        error('the compareAt input must be a scalar numeric value.')
    end
end

%% parse optional inputs
%ensure all remaining inputs are consistent with name-value pair format
if mod(length(varargin),2) > 0 || ...
        ~all(cellfun(@ischar,varargin(1:2:length(varargin))))
    error(['all inputs (after samples and parameterRequest) ' ...
        'must be submitted as name-value pairs.'])
elseif length(unique(varargin(1:2:length(varargin))))<length(varargin(1:2:length(varargin)))
    error('at least one property name (in the name-value pair inputs) is a duplicate.')
end

%set default values for the optional inputs
support = [];
xi = [];

varargin(1:2:length(varargin)) = ... %case insensitive (convert all to lower case)
    cellfun(@lower,varargin(1:2:length(varargin)),'uni',0);
for v = 1:2:length(varargin)
    switch varargin{v}
        %----------------------------------------------------------------%
        case 'support'
            support = varargin{v+1};
            %check the given support input
            if ischar(support)
                if ~ismember(support,{'unbounded','positive'})
                    error(['the support input must be either ''unbounded'', ' ...
                        '''positive'', or a two-element numeric vector.'])
                end
                if compareAt < 0 && strcmp(support,'positive')
                    error(['the comparison value (%.2f) and the support ' ...
                        '(''positive'') are not compatible.'],compareAt)
                end
            elseif isnumeric(support)
                if ~isvector(support) || ~isequal(length(support),2)
                    error(['the support input must be either ''unbounded'', ' ...
                        '''positive'', or a two-element numeric vector.']) 
                end
                if compareAt < support(1) || compareAt > support(2)
                    error(['the comparison value (%.2f) is outside the ' ...
                        'support interval bounds ([%.2f %.2f]).'], ...
                        compareAt,support(1),support(2))
                end
                tolerance = 0.005*diff(support);
                if compareAt < support(1)*(1+tolerance) || ...
                    compareAt > support(2)*(1-tolerance)
                    error(['the comparison value (%.2f) is within the support ' ...
                        'interval ([%.2f %.2f]), but is farily close to a bound. ' ...
                        'using this simple, generalized code is inadvisable. ' ...
                        'it is suggested to write your own script.'], ...
                        compareAt,support(1),support(2))
                end
            else
                error(['the support input must be either ''unbounded'', ' ...
                    '''positive'', or a two-element numeric vector.']) 
            end
        %----------------------------------------------------------------%
        case 'xi'
            xi = varargin{v+1};
            if ~isnumeric(xi) || ~isvector(xi) || ~all(diff(xi)>0)
                error('xi input must be a vector of increasing values.')
            end
        %----------------------------------------------------------------%
        otherwise
            error('the optional input name ''%s'' was not recognized.', ...
                varargin{v})
    end
end

%% check for MATLAB's kernel-smoothed density function
if ~exist('ksdensity','file')
    error(['this function requires the ksdensity function ' ...
        '(from MATLAB''s Statistics toolbox).'])
end

%% Savage-Dickey Bayes factor
%reformat chains
priorChains = priorChains(:);
posteriorChains = posteriorChains(:);

if isempty(xi) && isempty(support)
    ksOptions = {};
elseif ~isempty(xi) && isempty(support)
    ksOptions = {xi};
elseif isempty(xi) && ~isempty(support)
    ksOptions = {'support',support};
elseif ~isempty(xi) && ~isempty(support)
    ksOptions = {xi,'support',support};
end

%a smoothed posterior distribution
[f_post,x_post] = ksdensity(posteriorChains,ksOptions{:});
%a smoothed prior distribution
[f_prior,x_prior] = ksdensity(priorChains,ksOptions{:});

%find the value in each distribution closest to the comparison value
%... for the posterior
[~,ind_post] = min(abs(x_post - compareAt));
%... for the prior
[~,ind_prior] = min(abs(x_prior - compareAt));

%ratio of the densities at the comparison value is the bayes factor
density_post = f_post(ind_post);
density_prior = f_prior(ind_prior);
BF01 = density_post/density_prior;
% logBF01 = log(density_post) - log(density_prior);
BF10 = 1/BF01;

%% make plot
%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos(1:2) [480 360]*figScaling];
close(dumf.Number); %close dummy figure
%... and an axis
hold on

%underlay a vertical line at the comparison value
ymax = max(max(f_prior),max(f_post));
plot(compareAt*[1 1],ymax*2*[0 1],'color',0.85*[1 1 1])

%prior & posterior densities
priorLineOptions = {'linewidth',linePt*1.5,'color',[0 0 0],'linestyle','--'};
postLineOptions = {'linewidth',linePt*1.5,'color',[0 0 0],'linestyle','-'};
priorMarkerOptions = {'linewidth',linePt*1.5, ...
    'marker','o','markersize',markSz, ...
    'markeredgecolor','k','markerfacecolor',[1 1 1]};
postMarkerOptions = {'linewidth',linePt*1.5, ...
    'marker','o','markersize',markSz, ...
    'markeredgecolor','k','markerfacecolor',[0 0 0]};
plot(x_prior,f_prior,priorLineOptions{:})
plot(x_post,f_post,postLineOptions{:})
plot(compareAt,f_prior(ind_prior),priorMarkerOptions{:})
plot(compareAt,f_post(ind_post),postMarkerOptions{:})
%(dummies for legend)
h = NaN([2 1]);
h(1) = plot(x_prior(1),ymax*2,priorLineOptions{:},priorMarkerOptions{:});
h(2) = plot(x_post(1),ymax*2,postLineOptions{:},postMarkerOptions{:});
set(gca,'ylim',[0 ymax*1.1])

%format plot
ylabel('probability density')
set(gca,'box','on','fontsize',fontSz)
if BF10 > 1
    title(sprintf('BF_1_0 = %g',BF10), ...
        'interpreter','tex')
else
    title(sprintf('BF_1_0 = %.2g (BF_0_1 = %g)',BF10,BF01), ...
        'interpreter','tex')
end
legend(h,'prior','posterior','location','best','box','off')

end
