function MCMCtable = mcmctable(samples,method)
%MCMCTABLE collects diagnostic & posterior summary statistics in a table.
% 
% MCMCTABLE = MCMCTABLE(SAMPLES)
%   this function will create a table of statistics, MCMCTABLE.  
%   
%   for all instances of all parameters found in SAMPLES, selected summary
%   statistics and diagnostics are computed from the posterior samples,
%   collected in a table, and returned as SUMMARY. 
%   the SUMMARY table will contain the following variables:
%       'parameter'     >>  parameter names (scalar or indexed instance)
%       'mean'          >>  the mean of the posterior samples
%       'sd'            >>  the standard deviation of the posterior samples
%       'q50'           >>  the mean of the posterior samples
%       'mad'           >>  median absolute deviation of the posterior samples
%       'rhat           >>  folded split rhat
%       'ess_bulk'      >>  bulk effective sample size
%       'ess_tail'      >>  tail effective sample size
%   all diagnostic quantities are computed in accordance with current best
%   practice:   Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
%   as such, at time of writing (6/2021), this function's output matches
%   the summary output from Stan (v 2.26.1).
% 
% 
% MCMCTABLE = MCMCTABLE(SAMPLES,'BDA3')
%   instead of computing the currently accepted convergence diagnostics, 
%   one can elect to use the formulae of BDA3:
%       'rhat'      >>  split rhat
%     	'n_eff'     >>  effective sample size
%   these quantities are computed in accordance with the BDA3 text:
%               Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%                   Bayesian Data Analysis, 3rd ed. CRC.
%   this may be useful for backward compatability with previous work.
% 
% MCMCTABLE = MCMCTABLE(SAMPLES,PARAMETERREQUEST)
%   generates an MCMCTABLE output for the requested parameters only.
%   PARAMETERREQUEST must be a valid parameter name string or cell of 
%   valid parameter name strings. instance names may be included.
% 
% 
% See also INTERPRETDIAGNOSTICS, COMPUTERHAT, COMPUTEESS, EXTRACTSAMPLES
% 
% (c) beth baribault 2021 ---                                 > matstanlib 

%% check inputs
if nargin > 2
    error('too many inputs.')
end

%samples
if ~isstruct(samples)
    error(['first input must be a structure of posterior samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end

% %parameterRequest
% if nargin < 2 || isempty(parameterRequest)
%     parameterRequest = fieldnames(samples); %default
% else
%     if ischar(parameterRequest)
%         parameterRequest = {parameterRequest};
%     elseif ~iscell(parameterRequest) || ~all(cellfun(@ischar,parameterRequest))
%         error(['parameterRequest must be input as a string or ' ...
%                'cell of strings.'])
%     end
% end

%method
validMethods = {'current', ...      %folded-split Rhat, ESS_bulk, ESS_tail
                'BDA3', ...         %split Rhat, N_eff
                'BDA2', };          %Rhat
if nargin < 2
    method = 'current'; %default is current best practice
elseif ischar(method)
    if ~ismember(method,validMethods)
        error(['''%s'' is not a valid Rhat computation method.  \n' ...
            'valid computation methods include: ''%s'''], ...
            method,strjoin(validMethods,''', '''))
    end
else
    error('Rhat computation method must be a string (@ischar==true).')
end

%% calculate convergence diagnostics
%create a list of parameter instances
parameters = getparaminstances([], ...
    fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));
isInstance = cellfun(@(x) any(x=='['),parameters);
nParameters = length(parameters);

%%% [6/2021] cmdstan 2.26.1, stansummary
%%%'    Mean  MCSE  StdDev  5%  50%  95%  N_Eff  N_Eff/s  R_hat'

%%% [6/2021] posterior, summarise_draws(x)
%%% #>    variable  mean median    sd   mad      q5   q95  rhat ess_bulk ess_tail

%create a list of quantities to compute
switch method
    case 'current'
        varNames = {'mean','median','sd','mad','q5','q95', ...  %summary statistics
                    'rhat','ess_bulk','ess_tail'};              %convergence statistics
    case {'BDA3','BDA2'}
        varNames = {'mean','median','sd','mad','q5','q95', ...  %summary statistics
                    'rhat','n_eff'};                            %convergence statistics
end
nVariables = length(varNames);

%start a table
% MCMCtable = table('size',[nParameters nVariables], ...
%     'RowNames',parameters, ...
%     'VariableNames',varNames, ...
%     'VariableTypes',repelem({'double'},nVariables));
MCMCtable = array2table(NaN([nParameters nVariables]), ...
    'RowNames',parameters,'VariableNames',varNames);

for p = 1:length(parameters)
    %account for parameters vs. parameter instances
    if isInstance(p)
        [parameter,ind] = str2ind(parameters{p});
        chains = samples.(parameter)(:,:,ind{:});
    else
        parameter = parameters{p};
        chains = samples.(parameter);
    end
    
    %fill in summary statistics
    for q = 1:length(varNames)
        quantity = varNames{q};
        switch quantity
            case 'mean'
                value = mean(chains(:));
            case 'median'
                value = median(chains(:));
            case 'sd'
                value = std(chains(:),0); %normalize by N-1
            case 'mad'
                value = median(abs( chains(:) - median(chains(:)) ));
            case 'q5'
                value = quantile(chains(:),0.05);
            case 'q95'
                value = quantile(chains(:),0.95);
            case 'rhat'
                value = computerhat(chains,method);
            case 'ess_bulk'
                value = computeess(chains,'bulk');
            case 'ess_tail'
                value = computeess(chains,'tail');
            case 'n_eff'
                value = computeess(chains,'BDA3');
        end
        MCMCtable.(quantity)(p) = value;
    end
    clearvars chains 
end

end

%%% [6/2021]
%%% cmdstan version 2.26.1
%%% stansummary output still appears to use the BDA3 calculations
% 
% [stanSummaryTxt,stanSummary] = fit.print('sig_figs',5);
% mcmcSummary = mcmctable(samples,'BDA3');
% isequal(stanSummary.R_hat(8:end),round(mcmcSummary.rhat,5,'significant'))
% 
% ans =
% 
%   logical
% 
%    1