function essplot(samples,parameterRequest)
%ESSEFFICIENCYPLOT generates the "efficiency per iteration" diagnostic ESS plot.  
% 
% this function visualizes the "efficiency per iteration" plot recommended
% by Vehtari et al. (2020), one of three new diagnostic plots based on
% effective sample size (ESS) estimates.
% 
% if the marginal posterior distribution has been efficiently explored,
% then ESS estimates should grow linearly with the number of samples.  
% be wary if ESS estimates level off or decrease.  (some sampling
% inefficiencies may only become evident when sufficiently long chains are
% run.) 
% 
% the dashed line represents the minimum ESS, 100*(number of chains).
% 
% 
% ** NOTE: this function will generate no more than 25 figures per call. 
% 
% ESSEFFICIENCYPLOT(SAMPLES,PARAMETERREQUEST)
%   SAMPLES is assumed to be a strucutre of poserior samples in the format 
%   generated by extractsamples.m.  
%   
%   if PARAMETERREQUEST is a parameter name string, this function 
%   generates a figure with [all instances of] that parameter only.
%   if PARAMETERREQUEST is a cell of parameter name strings, this 
%   function generates a figure with [all instances of] each of those 
%   parameters.
%   PARAMETERREQUEST may include specific parameter instances, by 
%   including a valid index or indices in brackets after the parameter 
%   name (e.g., 'mu[3]','sigma[2,1]').  wildcards may also be used to 
%   select a subset of instances of a given parameter (e.g.,
%   'sigma[2,*]').
% 
% 
% Reference:  Vehtari, Gelman, Simpson, Carpenter, Bürkner (2020). 
%                   Rank-normalization, folding, and localization: An
%                   improved R^ for assessing convergence of MCMC.  ArXiv.
% 
% 
% See also COMPUTEESS, ESSPLOTS, PLOTMCSE, RANKPLOTS, TRACEDENSITY
% 
% (c) beth baribault 2021 ---                                 > matstanlib

msl.options

%% parse required inputs
if nargin < 2
    error('too few inputs.')
end

%samples
if ~isstruct(samples)
    error('first input must be a structure of posterior samples.')
else
    %extract size info immediately
    allFields = fieldnames(samples);
    [nIterations,nChains] = size(samples.(allFields{1}),[1 2]);
    nSamples = nIterations*nChains;
end

%parameterRequest
if isempty(parameterRequest)
    parameterRequest = fieldnames(samples);
elseif ~(iscell(parameterRequest) && all(cellfun(@ischar,parameterRequest)))
    if ischar(parameterRequest)
        parameterRequest = {parameterRequest};
    else
        error(['parameterRequest must be a string or cell of strings ' ...
               'containing valid parameter or instance name strings.'])
    end
end

%% prepare to generate plots
ESStypes = {'bulk','tail'};
ESScolors = getcolors('blue','lightblue');

%maximum number of figures to generate
maxNfigures = 25;

%create a list of parameter instances
parameters = getparaminstances(parameterRequest, ...
    fieldnames(samples),struct2cell(structfun(@size,samples,'uni',0)));

%warn if over the figure limit
if length(parameters) > maxNfigures
    warning(['requested %i traceplots, but this function will only ' ...
        'generate up to %i traceplots in a single call (to protect ' ...
        'against memory overload and subsequent freezing). ' ...
        'try requesting fewer parameters/parameter instances.  ' ...
        'continuing to generate the first %i traceplots now ...'], ...
        length(parameters),maxNfigures,maxNfigures)
end

%a quick dummy plot (to get nicer x-axis ticks)
dumf = figure;
plot(1:nSamples)
xt = get(gca,'xtick');
if xt(1)==0, xt(1) = []; end
close(dumf.Number)

%% make an "efficiency per iteration" plot for each parameter

nBatches   = 20; %maximum

%minimum quantites for plot
minNbatches = 4;
minSamplesPerBatch = 100;

%determine how many iterations to include in each batch
batchOK = false;
while ~batchOK
    nItersPerBatch = floor(nIterations/nBatches);
    if nItersPerBatch*nChains >= minSamplesPerBatch
        batchOK = true;
    else
        nBatches = nBatches - 2;
        if nBatches < minNbatches
            error('need at least %i samples total to make this plot.', ...
                minNbatches*minSamplesPerBatch)
        end
    end
end
batchIters = nIterations:-nItersPerBatch:0; %endpoints of equally-sized 
batchIters = batchIters(end:-1:1);          %batches ... in order
firstIter = batchIters(1)+1;        %iter to start ALL batches
lastIter =  batchIters(2:end);      %iter to end EACH batch
nBatches = length(lastIter);

%make ESS evolution plot
isInstance = cellfun(@(x) any(x=='['),parameters);
for p = 1:min([length(parameters),maxNfigures])
    %account for parameters vs. parameter instances
    if isInstance(p)
        [parameter,ind] = str2ind(parameters{p});
        chains = samples.(parameter)(:,:,ind{:});
    else
        parameter = parameters{p};
        chains = samples.(parameter);
    end
    
    %start a figure ...
    dumf = figure(999); %dummy figure to protect sizing
    f = figure('color',[1 1 1]);
    fpos = f.Position;
    f.Position = [fpos(1:2) [600 380]*figScaling];
    close(dumf.Number); %close dummy figure
    hold on
    
    %%% underlay lines at ESS criterion
    plot([1 nSamples],nChains*100*[1 1],'k--')
    
    %%% plot efficiency per iteration
    h = gobjects([1 length(ESStypes)]);
    ESSbyIter = NaN([1 nBatches]);
    for c = 1:length(ESStypes)
        ESStype = ESStypes{c};
        for n = 1:nBatches
            thisBatch = firstIter:lastIter(n);
            ESSbyIter(n) = computeess(chains(thisBatch,:),ESStype);
        end        
        h(c) = plot(lastIter*nChains,ESSbyIter,'-o', ...
            'color',ESScolors(c,:),'markerfacecolor',ESScolors(c,:), ...
            'linewidth',linePt*2);
    end
    
    %%% format
    xlim([1 nSamples])
    set(gca,'xtick',xt)
    yl = get(gca,'ylim'); yl = [0 max([yl(2) 1.1*100*nChains])]; set(gca,'ylim',yl)
    yt = get(gca,'ytick'); yt(yt<0) = []; set(gca,'ytick',yt)
    xlabel('total number of draws')
    ylabel('ESS')
    title(parameters{p},'interpreter','none')
    legend(h,ESStypes, ...
        'location','southoutside','orientation','horizontal')
    set(gca,'fontsize',fontSz,'box','on')
end

end