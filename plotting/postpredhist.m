function postpredhist(observedData,predictedData,parameter)
%POSTPREDHIST plots a unidimensional posterior predictive check histogram. 
% 
% POSTPREDHIST(OBSERVEDDATA,PREDICTEDDATA)
%   this function generates a single figure, in which smoothed posterior 
%   predictive data is overlaid on a histogram of the observed data. 
% 
% POSTPREDHIST(...,PARAMETER)
%   including PARAMETER, a parameter name string, allows for more 
%   descriptive axis labels, etc. 
% 
% See also PLOTDENSITY, SMOOTHDENSITY
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% parse inputs
if ismatrix(predictedData)
    nDatasets = 1;
else
    nDatasets = size(predictedData,3);
end
if nargin == 2
    parameter = 'y';
end

%% plot
%test discrete
% observedData = binornd(10,0.3,size(observedData));
% predictedData = binornd(10,0.35,size(predictedData));

%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) fpos([3 4])*figScaling];
close(dumf.Number); %close dummy figure
%... and an axes
ax = gca;
hold on
%... and some empty handles
h = gobjects([2 1]);

%plot the observed data as a histogram
obsIsDiscrete = all(~mod(observedData,1));
if obsIsDiscrete
    %bar plot
    [fObserved,xObserved] = smoothdensity(observedData);
    h(1) = bar(xObserved,fObserved,'barwidth',0.7,'linewidth',2, ...
        'facecolor',getcolors('lightgray'),'edgecolor',getcolors('gray'));
else
    h(1) = histogram('BinEdges',[0 1],'BinCounts',0, ...
        'facecolor',getcolors('lightgray'),'edgecolor',getcolors('gray'), ...
        'linewidth',2); %just for legend
    %histogram
    [fObserved,xObserved] = histcounts(observedData,'normalization','pdf');
    histogram('BinEdges',xObserved,'BinCounts',fObserved, ...
        'facecolor',getcolors('lightgray'),'edgecolor',getcolors('lightgray'));
    histogram('BinEdges',xObserved,'BinCounts',fObserved, ...
        'displaystyle','stairs','linewidth',2,'edgecolor',getcolors('gray'));
end

%plot the predicted data as a smoothed density
predMeanColor = getcolors('mutedblue');
predDatasetColor = mean([predMeanColor;1 1 1]);
predIsDiscrete = all(~mod(predictedData,1));
if predIsDiscrete
    xPredicted = unique(predictedData(:))';
    fPredicted = NaN([nDatasets length(xPredicted)]);
    for n = 1:nDatasets
        predictedDataSet = predictedData(:,:,n);
        f = histcounts(predictedDataSet, ...
            [xPredicted-0.5 xPredicted(end)+0.5], ...
            'normalization','pdf');
        fPredicted(n,:) = f;
    end
    %circular markers
    for n = 1:nDatasets
        plot(xPredicted,fPredicted,'linestyle','none', ...
            'marker','o','markersize',markSz*1,5, ...
            'markerfacecolor','none','markeredgecolor',predDatasetColor);
    end
    fPredicted = mean(fPredicted,1);
    h(2) = plot(xPredicted,fPredicted,'linestyle','none', ...
        'marker','o','markersize',markSz*1.5, ...
        'markerfacecolor',predMeanColor,'markeredgecolor',predMeanColor);
else %predIsContinuous
    xPredicted = linspace(min(predictedData(:)),max(predictedData(:)), ...
        min(size(predictedData,1),1000));
    fPredicted = NaN([nDatasets length(xPredicted)]);
    for n = 1:nDatasets
        predictedDataSet = predictedData(:,:,n);
        f = ksdensity(predictedDataSet(:),xPredicted);
        fPredicted(n,:) = f;
    end
    %smooth line
    for n = 1:nDatasets
        plot(xPredicted,fPredicted(n,:),'linewidth',linePt*0.75,'color',predDatasetColor)
    end
    fPredicted = mean(fPredicted,1);
    h(2) = plot(xPredicted,fPredicted,'linewidth',linePt*1.25,'color',predMeanColor);
    ax.XLim = [min(xPredicted) max(xPredicted)];
end

%format plot
xLimits = [min(min(xObserved),min(xPredicted)) ...
           max(max(xObserved),max(xPredicted))];
ax.XLim = [xLimits(1) - 0.1*diff(xLimits) ...
           xLimits(2) + 0.1*diff(xLimits)];
ax.YLim = [0 max(max(fObserved),max(fPredicted))*1.1];
xlabel([parameter '~'])
ylabel(sprintf('p(%s~|%s)',parameter,parameter))
legend(h,sprintf(' %s: observed data',parameter), ...
         sprintf('%s~: predicted data',parameter))
set(ax,'box','on')

end
