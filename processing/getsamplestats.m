function estimatedValues = getsamplestats(samples,varargin)
%GETSAMPLESTATS computes common summary statistics for posterior samples. 
% 
% ESTIMATEDVALUES = GETSAMPLESTATS(SAMPLES)
%   this function returns a structure (ESTIMATEDVALUES) with one fieldname
%   for each parameter (i.e., for each fieldname found in SAMPLES).  
%   for each instance of each parameter, commonly-used summary statistics
%   will be computed from the relevant posterior samples.
%   
%   ESTIMATEDVALUES has the following structure:
%       estimatedValues 
%           .[parameter]  >> parameter name (e.g., beta)
%               .name       >> instance name strings (e.g., 'beta[56]')
%               .mean       >> mean
%               .median     >> median
%               .mode       >> mode
%               .lowerCI95  >> 95% central credible interval lower bound
%               .upperCI95  >> 95% central credible interval upper bound
%               .lowerCI90  >> 90% central CI lower bound
%               .upperCI90  >> 90% central CI upper bound
%               .lowerCI50  >> 50% central CI lower bound
%               .upperCI50  >> 50% central CI upper bound
%   note that each subfield's value (e.g., estimatedValues.beta.mean) will
%   be of the same size as the relevant parameter. 
% 
% 
% the following optional inputs may be given in any order as they are
% recognized by type: 
% 
% ESTIMATEDVALUES = GETSAMPLESTATS(...,TRUEVALUES)
%   TRUEVALUES is a structure containing the actual, known parameter values
%   used to simulate the data, with one field per parameter.  
%   if TRUEVALUES is given, then the ESTIMATEDVALUES structure will contain
%   two additional subfields:
%       estimatedValues
%           .[parameter]
%               .true   >> the true parameter value, used to simulate
%                               the data 
%               .inCI95 >> whether the true parameter value is within the 
%                           95% central CI (1 if so; 0 if not)
%               .inCI90 >> whether the true parameter value is within the 
%                           90% central CI (1 if so; 0 if not)
%               .inCI50 >> whether the true parameter value is within the 
%                           50% central CI (1 if so; 0 if not)
% 
% ESTIMATEDVALUES = GETSAMPLESTATS(...,'HDI')
%   if the string 'HDI' or 'HPD' is included in the inputs, then various
%   highest posterior density intervals will be computed also.
%   multiple addtional subfields will be added to ESTIMATEDVALUES: 
%       estimatedValues
%           .[parameter]
%               .lowerHDI95 >> 95% highest density interval lower bound
%               .upperHDI95 >> 95% highest density interval upper bound
%               .lowerHDI90 >> 90% HDI lower bound
%               .upperHDI90 >> 90% highest density interval upper bound
%               .lowerHDI50 >> 50% highest density interval lower bound
%               .upperHDI50 >> 50% highest density interval upper bound
%               .inHDI95    >> [if applicable] whether the true parameter
%                               value is within the 95% HDI
%               .inHDI90    >> [if applicable] whether the true parameter
%                               value is within the 90% HDI
%               .inHDI50    >> [if applicable] whether the true parameter
%                               value is within the 50% HDI
% 
% ESTIMATEDVALUES = GETSAMPLESTATS(...,PARAMETERREQUEST)
%   computes statistics for a given parameter or parameters only. 
%   parameter names must be exact.
% 
% NOTE: this function will calculate statistics for all instances of 
% each [requested] parameter. the dimension of the calculated statistics 
% will match the parameter's dimensions, such that for a parameter 
% whose field in SAMPLES is of size:
%     [nIterations nChains <arrayDimensions> <parameterDimensions>]
% each field in ESTIMATEDVALUES will be of size:
%     [<arrayDimensions> <parameterDimensions>] 
% for example, if the parameter 'mu' is a 3-element vector, then
% each field in estimatedValues.mu will contain a 3-element vector. 
% 
% 
% See also EXTRACTSAMPLES, PLOTRECOVERY, COMPUTECREDINT, COMPUTEHDI
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% check inputs
if ~isstruct(samples)
    error(['first input must be a structure of posterior samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end

trueValuesKnown = false;
includeHDI = false;
for v = 1:length(varargin)
    %true values structure
    if isstruct(varargin{v})
        t = varargin{v};
        trueValuesKnown = true;
    %include HDI
    elseif ischar(varargin{v})
        if strcmpi('hdi',varargin{v}) || strcmpi('hpd',varargin{v}) %case insensitive
            includeHDI = true;
        else
    %parameter request
            parameterRequest = {varargin{v}};
        end
    elseif iscell(varargin{v}) && all(cellfun(@ischar,varargin{v}))
        parameterRequest = varargin{v};
    else
        error(['input type not recognized. ' ...
            'parameter names must be input as a string or cell of ' ...
            'strings. true values must be input as a struct. ' ...
            'no other inputs are accepted.'])
    end
end
if nargin > 5
    error('invalid number of input arguments.')
end

%% extract [requested] parameter names
parameters = fieldnames(samples);
if exist('parameterRequest','var')
    if any(~ismember(parameterRequest,parameters))
        error(['the parameter request includes at least one string ' ...
            'that is not a valid parameter name.'])
    end
    parameters = parameters(ismember(parameters,parameterRequest));
end

%% calculate stats
for p = 1:length(parameters)
    parameter = parameters{p};
    chainSize = size(samples.(parameter));
    if length(chainSize) == 2
        %scalar parameter
        parameterSize = [1 1];
    elseif length(chainSize) == 3
        %vector parameter
        parameterSize = [chainSize(end) 1];
    else
        %N-dimensional matrix parameter
        parameterSize = chainSize(3:end);
    end
    instanceNames = getparaminstances(parameter, ...
        fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));
    estimatedValues.(parameter).name = instanceNames;
    %preallocate
    estimatedValues.(parameter).true    = NaN;
    estimatedValues.(parameter).mean    = NaN(parameterSize);
    estimatedValues.(parameter).median  = NaN(parameterSize);
    estimatedValues.(parameter).mode    = NaN(parameterSize);
    estimatedValues.(parameter).argmax  = NaN(parameterSize);
    estimatedValues.(parameter).lowerCI95 = NaN(parameterSize);
    estimatedValues.(parameter).upperCI95 = NaN(parameterSize);
    estimatedValues.(parameter).lowerCI90 = NaN(parameterSize);
    estimatedValues.(parameter).upperCI90 = NaN(parameterSize);
    estimatedValues.(parameter).lowerCI50 = NaN(parameterSize);
    estimatedValues.(parameter).upperCI50 = NaN(parameterSize);
    if includeHDI
        estimatedValues.(parameter).lowerHDI95 = NaN(parameterSize);
        estimatedValues.(parameter).upperHDI95 = NaN(parameterSize);
        estimatedValues.(parameter).lowerHDI90 = NaN(parameterSize);
        estimatedValues.(parameter).upperHDI90 = NaN(parameterSize);
        estimatedValues.(parameter).lowerHDI50 = NaN(parameterSize);
        estimatedValues.(parameter).upperHDI50 = NaN(parameterSize);
    end
    if trueValuesKnown
        thisParamKnown = false;
        if isfield(t,parameter)
            trueValue = t.(parameter);
            if isequal(size(trueValue),parameterSize)
                estimatedValues.(parameter).true = trueValue;
                estimatedValues.(parameter).inCI95 = NaN(parameterSize);
                estimatedValues.(parameter).inCI90 = NaN(parameterSize);
                estimatedValues.(parameter).inCI50 = NaN(parameterSize);
                if includeHDI
                    estimatedValues.(parameter).inHDI95 = NaN(parameterSize);
                    estimatedValues.(parameter).inHDI90 = NaN(parameterSize);
                    estimatedValues.(parameter).inHDI50 = NaN(parameterSize);
                end
                thisParamKnown = true;
            else
%                 warning(['skipping the ''true'' and ''inCI/inHDI'' subfields ' ...
%                     'for %s (due to a mismatch between the size of %s and ' ...
%                     'the size of the true values for %s).'], ...
%                     parameter,parameter,parameter)
            end
        else
%             warning(['skipping the additional subfields for %s.' ...
%                 '(because while there are samples for %s, %s was ' ...
%                 'not included in the true values structure).'], ...
%                 parameter,parameter,parameter)
        end
    end
    %calculate posterior statistics
    nInstances = length(instanceNames);
    clearvars ind
    % *** loop over parameters
    for n = 1:nInstances
        [ind{1:length(parameterSize)}] = ind2sub(parameterSize,n);
        if nInstances == 1
            chains = samples.(parameter);
        else
            chains = samples.(parameter)(:,:,ind{:});
        end
        chains = chains(:);
        estimatedValues.(parameter).mean(ind{:})    = mean(chains);
        estimatedValues.(parameter).median(ind{:})  = median(chains);
        estimatedValues.(parameter).mode(ind{:})    = mode(chains);
        try 
            [f,x] = smoothdensity(chains);
            estimatedValues.(parameter).argmax(ind{:}) = x(f==max(f));
        catch
            %skip
        end
        estimatedValues.(parameter).lowerCI95(ind{:}) = quantile(chains,0.025);
        estimatedValues.(parameter).upperCI95(ind{:}) = quantile(chains,0.975);
        estimatedValues.(parameter).lowerCI90(ind{:}) = quantile(chains,0.05);
        estimatedValues.(parameter).upperCI90(ind{:}) = quantile(chains,0.95);
        estimatedValues.(parameter).lowerCI50(ind{:}) = quantile(chains,0.25);
        estimatedValues.(parameter).upperCI50(ind{:}) = quantile(chains,0.75);
        if includeHDI
            hdi = computecredint(chains,0.95,'hdi');
            estimatedValues.(parameter).lowerHDI95(ind{:}) = hdi(1);
            estimatedValues.(parameter).upperHDI95(ind{:}) = hdi(2);
            hdi = computecredint(chains,0.90,'hdi');
            estimatedValues.(parameter).lowerHDI90(ind{:}) = hdi(1);
            estimatedValues.(parameter).upperHDI90(ind{:}) = hdi(2);
            hdi = computecredint(chains,0.50,'hdi');
            estimatedValues.(parameter).lowerHDI50(ind{:}) = hdi(1);
            estimatedValues.(parameter).upperHDI50(ind{:}) = hdi(2);
        end
        if trueValuesKnown && thisParamKnown
            estimatedValues.(parameter).inCI95(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerCI95(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperCI95(ind{:});
            estimatedValues.(parameter).inCI90(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerCI90(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperCI90(ind{:});
            estimatedValues.(parameter).inCI50(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerCI50(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperCI50(ind{:});
            if includeHDI
              estimatedValues.(parameter).inHDI95(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerHDI95(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperHDI95(ind{:});
              estimatedValues.(parameter).inHDI90(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerHDI90(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperHDI90(ind{:});
              estimatedValues.(parameter).inHDI50(ind{:}) = ...
                trueValue(ind{:}) >= estimatedValues.(parameter).lowerHDI50(ind{:}) & ...
                trueValue(ind{:}) <= estimatedValues.(parameter).upperHDI50(ind{:});
            end
        end
    end
end
    
end