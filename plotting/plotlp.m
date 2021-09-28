function plotlp(diagnostics,rtable)
%PLOTLP summarizes some diagnostic metrics in a single plot.
% 
% PLOTLP(DIAGNOSTICS,RTABLE)
%   generates a figure with three subplots based on a DIAGNOSTICS, a 
%   structure of Stan-generated diagnostic quantities, and RTABLE, a 
%   table of calculated diagnostic quantities. 
%   on the left of the figure is a trace of the total log posterior 
%   density, with a rug plot of divergent transitions overlaid at the 
%   bottom.  
%   on the right of the figure is a bar plot summary of rhat values.  
% 
% See also PLOTAUTOCORR, TRACEDIVERGENT, RHATTABLE, EXTRACTSAMPLES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

msl.options

%% parse required inputs
if ~isstruct(diagnostics) || ~isfield(diagnostics,'divergent__')
    error(['first input must be a structure of diagnostic quantities ' ...
        '(consisitent with the output of extractsamples.m).'])
end
if ~istable(rtable) || ~ismember('rhat',rtable.Properties.VariableNames)
    error(['second input must be a table of rhats and other values' ...
        '(consisitent with the output of rhattable.m).'])
end

%%
%extract log posterior probability
lpChains = diagnostics.lp__;
%extract divergent transitions
divergences = diagnostics.divergent__;
%extract dimensions
nIterations = size(lpChains,1);
nChains = size(lpChains,2);

%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos(1:2) [800 420]*figScaling];
close(dumf.Number); %close dummy figure
%... and an axis
ax1 = axes('position',[0.09  0.35  0.55  0.60]);
ax2 = axes('position',[0.09  0.09  0.55  0.20]);
ax3 = axes('position',[0.74  0.275 0.235 0.65]);

%get some colors
chainColor = getcolors('mcmc');

%% plot lp ...
axes(ax1)
hold on
h = gobjects([1 nChains]);
for n = 1:nChains
    plot(lpChains(:,n),'color',chainColor(n,:),'linewidth',linePt)
end
%dummy lines for legend
for n = 1:nChains, h(n) = plot(-1,lpChains(1,1),'color',chainColor(n,:),'linewidth',4); end
legend(h,arrayfun(@num2str,1:nChains,'uni',0), ...
    'location','southoutside','orientation','horizontal')
%format plot
set(gca,'xlim',[0 nIterations+1])
set(gca,'fontsize',fontSz,'box','on')
ylabel('total log probability density','fontsize',fontSz)

%% ... and add divergences as a rugplot at the bottom
axes(ax2)
hold on, axis ij
leftLabels = cell([nChains 1]);
rightLabels = cell([nChains 1]);
for n = 1:nChains
    divSamples = find(divergences(:,n)==1);
    if ~isempty(divSamples)
        for c = 1:length(divSamples)
            plot(divSamples(c)*[1 1],n+[-0.4 0.4], ...
                'linewidth',linePt*2,'color',getcolors('red'))
        end
    end
    leftLabels{n} = sprintf('chain %i',n);
    rightLabels{n} = num2str(length(divSamples));
end
%format
set(gca,'xticklabels',{})
xlim([0 nIterations + 1])
ylim([0 nChains + 1])
%left y axis
yyaxis left
set(gca,'ylim',[0 nChains + 1],'ytick',1:nChains,'yticklabels',leftLabels)
set(gca,'ygrid','on','ycolor','k')
%right y axis
yyaxis right, axis ij
set(gca,'ylim',[0.25 nChains + 0.75],'ytick',1:nChains,'yticklabels',rightLabels)
set(gca,'ycolor','k')
set(gca,'fontsize',fontSz,'box','on')
xlabel('divergent transitions','fontsize',fontSz)
set(gca,'xcolor','w')

%% bar plot of rhat
axes(ax3)
hold on
%bins for rhat
rhatEdges = [0,1.01,1.05,1.1,1.5,5,Inf];
nBins = length(rhatEdges) - 1;
rhatCounts = histcounts(rtable.rhat,rhatEdges);

%get a nice color map
nGoodBins = sum(rhatEdges < 1.1);
m1 = makecolormap(getcolors('darkgray'),getcolors('gray'),nGoodBins);
m2 = makecolormap(getcolors('gray'),getcolors('red'),nBins-nGoodBins+1);
barColors = [m1;m2(2:end,:)];
%bar plot
h = gobjects([nBins 1]);
rhatCounts(rhatCounts==0) = eps;
textpos = max(rhatCounts)*0.1;
for b = 1:nBins
    h(b) = bar(b,rhatCounts(b),'facecolor',barColors(b,:));
    text(b,rhatCounts(b)+textpos,sprintf('%i',round(rhatCounts(b))), ...
        'horizontalalignment','center','fontsize',fontSz-2)
end
set(gca,'xlim',[0.5,nBins+0.5], ...
    'xtick',1:nBins,'xticklabel',{},'tickdir','out')
set(gca,'ylim',[0 max(rhatCounts)+textpos*2])
catEdges =  arrayfun(@num2str,rhatEdges(2:end-1),'uni',0);
categories = cell([nBins 1]);
for c = 2:length(categories)-1
    categories{c} = [catEdges{c-1} ' to ' catEdges{c}];
end
categories{1} = ['< ' catEdges{1}];
categories{end} = ['>' catEdges{end}];
set(gca,'xticklabel',categories,'xticklabelrotation',45)
set(gca,'fontsize',fontSz,'box','off')
xlabel('$\hat{R}$','fontsize',fontSz+4,'interpreter','latex')
ylabel('number of parameters')

end