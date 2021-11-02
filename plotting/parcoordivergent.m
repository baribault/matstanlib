function parcoordivergent(samples,diagnostics,parameterRequest)
%PARCOORDIVERGENT creates a parallel coordinate plot and highlights divergences. 
% 
% this function generates a single plot in which each line represents the 
% values of all parameter values for a single mcmc sample.  
% values sampled during divergent transitions are highlighted in red.  
% this plot may be useful in assessing the cause of divergent transitions 
% (and therefore altering a model such that it no longer generates 
% divergences).
% 
% PARCOORDIVERGENT(SAMPLES,DIAGNOSTICS)
%   SAMPLES is a structure containing mcmc samples and DIAGNOSTICS is a
%   structure containing Stan-generated diagnostics (both consistent 
%   with the format generated by extractsamples.m).
%   SAMPLES and DIAGNOSTICS are both requrired inputs.
% 
% NOTE: a maximum of 20 parameter instances may be included in a single
% plot. if parameters of interest are on very different scales, try 
% plotting pairs of parameters instead of more numerous groups. 
% 
% PARCOORDIVERGENT(SAMPLES,DIAGNOSTICS,PARMAETERREQUEST)
%   it will usually be useful (& necessary) to call this function with the
%   optional input PARAMETERREQUEST, a cell of parameter instance name 
%   strings.  if a parameter name string is included in PARAMETERREQUEST,
%   it will be expanded to parameter instance names.  
% 
% See also PLOTLP, TRACEDIVERGENT, EXTRACTSAMPLES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

msl.options

%% parse inputs
if ~isstruct(samples)
    error(['first input must be a structure of MCMC samples ' ...
        '(consisitent with the output of extractsamples.m).'])
end
if ~isstruct(diagnostics)
    error(['second input must be a structure of diagnostic quantities ' ...
        '(consisitent with the output of extractsamples.m).'])
end
if nargin == 2
    parameterRequest = fieldnames(samples);
elseif nargin > 2
    if isempty(parameterRequest)
        parameterRequest = fieldnames(samples);
    elseif ~(iscell(parameterRequest) && all(cellfun(@ischar,parameterRequest)))
        if ischar(parameterRequest)
            parameterRequest = {parameterRequest};
        else
            error(['parameter names must be input as a string or ' ...
                   'cell of strings.'])
        end
    end
end

%% create a list of parameter instances
parameters = getparaminstances(parameterRequest, ...
    fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));
nParameters = length(parameters);
isInstance = cellfun(@(x) any(x=='['),parameters);

%don't plot more than 20 parameters
if length(parameters) == 1
    error('too few parameters to create a plot.')
elseif length(parameters) > 20
    error(['too many parameters to plot (> 20).  try changing the ' ...
        'parameterRequest input to be more restrictive.'])
end

%% make parallel coordinates plot
%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) [100+60*nParameters 400]*figScaling];
close(dumf.Number); %close dummy figure
%... and an axis
hold on

%plot parameter values with divergent transitions highlighted
for nplot = 1:2
    if nplot == 1
        %first, plot good samples
        [iter,chain] = find(~diagnostics.divergent__);
        color = getcolors('black');
    elseif nplot == 2
        %next, plot divergent samples
        [iter,chain] = find(diagnostics.divergent__);
        color = getcolors('red');
    end
    for n = 1:length(iter)
        %select current sample
        sample = {iter(n),chain(n)};
        %get all parameter values for the current sample
        valuesAtSample = NaN([nParameters 1]);
        for p = 1:nParameters
            %account for parameters vs. parameter instances
            if isInstance(p)
                [parameter,ind] = str2ind(parameters{p});
                valuesAtSample(p) = samples.(parameter)(sample{:},ind{:});
            else
                parameter = parameters{p};
                valuesAtSample(p) = samples.(parameter)(sample{:});
            end
        end
        %plot all parameter values at current sample
        plot(1:nParameters,valuesAtSample,'color',color)
    end
end
%format plot
set(gca,'xlim',[0.9 nParameters+0.1], ...
    'xtick',1:nParameters,'xticklabel',parameters, ...
    'ticklabelinterpreter','none')
set(gca,'ygrid','on','box','on','fontsize',fontSz)

end