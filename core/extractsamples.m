function [samples,diagnostics] = extractsamples(interface,fit,varargin)
%EXTRACTSAMPLES extracts and reformats output from various Stan interfaces.
% 
% this function will extract and reformat the output from a variety of
% interfaces to Stan, such that the output is compatible with matstanlib.  
% 
% [SAMPLES,DIAGNOSTICS] = EXTRACTSAMPLES(INTERFACE,...)
%   the exact inputs required depend on which INTERFACE to Stan was used. 
%       'MATLABStan'  >>  https://github.com/brian-lau/MatlabStan
%       'Trinity'     >>  https://github.com/joachimvandekerckhove/trinity
%   INTERFACE is case insensitive.
%   
%   posterior samples are returned in the SAMPLES structure.  
%   diagnostic quantities are returned in the DIAGNOSTICS structure.  
% 
% [SAMPLES,DIAGNOSTICS] = EXTRACTSAMPLES('matlabstan',FIT)
%   if you have used MATLABStan as your interface to Stan, the StanFit
%   object, FIT, should be given.
% 
% [SAMPLES,DIAGNOSTICS] = EXTRACTSAMPLES('trinity',CHAINS,INFO)
%   if you have used Trinity as your interface to Stan, the CHAINS
%   structure and the INFO structure should both be given. 
%   posterior samples for all monitored parameters will be returned.
%   
%   regardless of INTERFACE, the posterior samples and diagnostic
%   quantities will be extracted, reformatted, and returned in a format 
%   that is compatible with all functions in the matstanlib library.  
%   
%   SAMPLES will have one field per parameter, and each field will contain
%   a matrix of posterior samples of size:  
%       [nIterations nChains <arrayDimensions> <parameterDimensions>]
%   (this is consistent with other common interfaces to Stan (e.g., PyStan)
%   but dissimilar to MATLAB interfaces for other Bayesian sampling engines
%   (such as MatJags, Trinity, etc.), in that instances of the same
%   parameter (e.g., mu[1] and mu[2]) are NOT assigned different fields in
%   the output structure.) 
%   
%   diagnostic quantities are returned in the DIAGNOSTICS structure.  
%   DIAGNOSTICS will have one field for each of Stan's monitored
%   diagnostics.  the fieldname for each quanitity includes the __ suffix
%   (e.g., divergent__).   
%   each field's value will be a [nIterations nChains]-sized matrix.  
% 
% 
% [SAMPLES,DIAGNOSTICS] = EXTRACTSAMPLES(...,PARAMETERREQUEST)
%   optionally, only the posterior samples for a requested set of
%   parameters, PARAMETERREQUEST, may be extracted and returned as SAMPLES. 
%   PARAMETERREQUEST must be a cell of valid parameter name strings.
%   parameter instance names are not accomodated. 
% 
% 
% See also GETSAMPLESTATS, RHATTABLE, REMOVECHAIN, MANUALBURN, 
%   EXPANDPARAMNAMES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 


%% check inputs
%interface
validInterfaces = {'matlabstan','trinity'};
if nargin < 1 || ~ischar(interface)
    error(['the first input must be a valid interface name string. ' ...
        'valid Stan interface names include: ''%s'''], ...
        strjoin(validInterfaces,''', '''))
else
    interface = lower(interface);
    if ~ismember(interface,validInterfaces)
        error(['the interface input ''%s'' was not recognized.  the only ' ...
            'valid Stan interface names include: ''%s'''], ...
            interface,strjoin(validInterfaces,''', '''))
    end
end

%check required inputs, given interface
switch interface
    % ------------------------------------------------------------------- %
    case 'matlabstan'
        if nargin < 2
            error(['a StanFit object is required to process output from ' ...
                'the MATLABStan interface to Stan.  your call to this ' ...
                'function should look like:%s\n' ...
                '    extractsamples(''matlabstan'',fit)'],'')
        end
    % ------------------------------------------------------------------- %
    case 'trinity'
        if nargin < 3
            error(['both a ''chains'' structure and an ''info'' structure ' ...
                'are required to process output from the Trinity interface. ' ...
                'your call to this function should look like:\n' ...
                '    extractsamples(''trinity'',chains,info)'])
        end
        %going to leave chains structure as 'fit' to avoid memory waste
        if ~isstruct(fit) || ~isstruct(varargin{1})
            error(['both a ''chains'' structure and an ''info'' structure ' ...
                'are required to process output from the Trinity interface, ' ...
                'but the second and/or third input is not a structure.'])
        end
        info = varargin{1};
        varargin = varargin(2:end);
        try
            info.options.monitorparams;
            info.options.allowunderscores;
        catch
            error(['the third input should be Trinity''s ''info'' ' ...
                'structure, but it doesn''t have some of the fields it ' ...
                'should have.  please check your inputs and try again.'])
        end
        if isempty(info.options.monitorparams)
            error(['''monitorparams'' must used as in input to callbayes' ...
                'in order for Trinity output to be compatible with matstanlib ' ...
                '(otherwise parameter names, etc. cannot be recovered).'])
        end
        if info.options.allowunderscores
            error(['callbayes'' ''allowunderscores'' option must be false ' ...
                'in order for Trinity output to be compatible with matstanlib ' ...
                '(otherwise parameter names, etc. cannot be recovered).'])
        end
end

%%% optional inputs %%%
%parameterRequest
if isempty(varargin)
    parameterRequest = {};
else
    parameterRequest = varargin{1};
    if ischar(parameterRequest)
        parameterRequest = {parameterRequest};
    elseif ~iscell(parameterRequest) || ~all(cellfun(@ischar,parameterRequest))
        error('parameterRequest must be a cell of parameter name strings.')
    end
end

%% extract & reformat samples
switch interface
  % --------------------------------------------------------------------- %
  case 'matlabstan'
    %% find all monitored quantites
    monitored = fieldnames(fit.sim.samples);
    %determine which are parameters vs. MCMC diagnostic metrics
    monitored_min2chars = cellfun(@(x) ['_' x],monitored,'uni',0);
    isParameter = ~cellfun(@(x) strcmp(x(end-1:end),'__'),monitored_min2chars);
    parameters = monitored(isParameter);
    
    %% accomodate parameter request
    if ~isempty(parameterRequest)
        if any(~ismember(parameterRequest,parameters))
            warning(['posterior samples for the following parameters were ' ...
                'requested, but these parameters were not monitored ' ...
                '(i.e., no samples for these parameters were found in ' ...
                'the MATLABStan output): ''%s'''], ...
                strjoin(parameterRequest(~ismember(parameterRequest,parameters))))
            error('invalid entries in parameterRequest input.')
        end
        parameters = parameters(ismember(parameters,parameterRequest));
    end
    
    %% extract samples
    %get core output dimensions
    nChains = length(fit.sim.samples);
    nIterations = size(fit.sim.samples(1).(parameters{1}),1);
    
    %extract samples for all monitored parameters
    for p = 1:length(parameters)
        %get dimensions of the current parameter
        dims = size(fit.sim.samples(1).(parameters{p}));
        dims = dims(2:end);
        %extracting all chains for the current parameter
        if dims == 1
            %if scalar, ...
            samples.(parameters{p}) = [fit.sim.samples.(parameters{p})];
        else
            %otherwise, ...
            %NOTE: array indices come BEFORE parameter indices!
            argDimsY = num2cell(dims);
            tmp = NaN(nIterations,nChains,argDimsY{:});
            colonDimsY(1:length(dims)) = {':'};
            for ch = 1:nChains
                tmp(:,ch,colonDimsY{:}) = fit.sim.samples(ch).(parameters{p});
            end
            samples.(parameters{p}) = tmp;
        end
    end
    
    %extract samples for all monitored diagnostic metrics
    metrics = monitored(~isParameter);
    for m = 1:length(metrics)
        diagnostics.(metrics{m}) = [fit.sim.samples.(metrics{m})];
    end
    
    
  % --------------------------------------------------------------------- %
  case 'trinity'
    parameters = info.options.monitorparams;
    % [~,sortOrder] = sort(cellfun(@length,parameters),'descend');
    % parameters = parameters(sortOrder)
    
    %% accomodate parameter request
    if ~isempty(parameterRequest)
        if any(~ismember(parameterRequest,parameters))
            warning(['posterior samples for the following parameters were ' ...
                'requested, but these parameters were not monitored ' ...
                '(i.e., no samples for these parameters were found in ' ...
                'the Trinity output): ''%s'''], ...
                strjoin(parameterRequest(~ismember(parameterRequest,parameters))))
        end
        parameters = parameters(ismember(parameters,parameterRequest));
    end
    
    %% recover parameter names and dimensions
    %get trinity instances --- with bad indexing
    trinityInsts = fieldnames(fit);
    
    %get dimensions
    nIterations = size(fit.(trinityInsts{1}),1);
    nChains = size(fit.(trinityInsts{1}),2);
    
    %for each parameter, starting with the longest name
    paramSize = cell(size(parameters));
    paramInstances = cell(size(parameters));
    paramIndices = cell(size(parameters));
    for p = 1:length(parameters)
        param = parameters{p};
        if sum(ismember(trinityInsts,param)) == 1
            %scalar parameter
            isScalarParam = cellfun(@(x) strcmp(x,param),trinityInsts);
            paramSize{p} = 0;
            paramInstances{p} = isScalarParam;
        else
            %multiple instances of this parameter
            isInstanceOfParam = cellfun(@(x) strncmp(x,[param '_'],length(param)+1),trinityInsts);
            if ~any(isInstanceOfParam)
                %parameter not found!
                paramSize{p} = {};
                warning(['Trinity''s ''info'' structure suggests that ' ...
                    'parameter ''%s'' was requested to be monitored, ' ...
                    'monitored, but samples for this parameter could ' ...
                    'not be found in the Trinity output.'],param)
            else
                %inferring size of this parameter
                instances = ... %add an underscore to all instances of this param
                    cellfun(@(x) [x '_'], trinityInsts(isInstanceOfParam),'uni',0);
                %extract indices from each instance
                nInstances = length(instances);
                nIndices = length(find(instances{1}=='_')) - 1;
                paramIndices{p} = NaN([nInstances nIndices]);
                for row = 1:nInstances
                    underscores = find(instances{row}=='_');
                    for col = 1:nIndices
                        paramIndices{p}(row,col) = str2double( ...
                            instances{row}(underscores(col)+1:underscores(col+1)-1));
                    end
                end
                paramSize{p} = max(paramIndices{p});
                paramInstances{p} = isInstanceOfParam;
            end
        end
    end    
    
    notFound = cellfun(@iscell,paramSize);
    parameters(notFound) = [];
    paramSize(notFound) = [];
    paramInstances(notFound) = [];
    paramIndices(notFound) = [];
    
    %% extract samples
    %extract samples for all monitored parameters
    for p = 1:length(parameters)
        if paramSize{p} == 0
            %scalar parameter
            samples.(parameters{p}) = fit.(parameters{p});
        else
            %multiple instances of this parameter
            samples.(parameters{p}) = NaN([nIterations nChains paramSize{p}]);
            instances = trinityInsts(paramInstances{p});
            for n = 1:length(instances)
                ind = num2cell(paramIndices{p}(n,:));
                samples.(parameters{p})(:,:,ind{:}) = fit.(instances{n});
            end
        end
    end
    
    %extract samples for all monitored diagnostic metrics
    diagnostics = info.tuning.chain;
end

end
