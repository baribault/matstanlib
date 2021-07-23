function plotdivergences(diagnostics)
%PLOTDIVERGENCES creates a rug plot of divergent transitions chain-by-chain.
% 
% PLOTDIVERGENCES(DIAGNOSTICS)
%   generates a figure with a rugplot indicating whether each sample in
%   each chain was the result of a divergent transition (red mark).
% 
% PLOTDIVERGENCES(DIVERGENCES)
%   alternatively, just the [nIterations nChains] matrix found in
%   diagnostics.divergent__ may be input.
% 
% See also PLOTAUTOCORR, TRACEDIVERGENT, RHATTABLE, EXTRACTSAMPLES
% 
% (c) beth baribault 2020 ---                                 > matstanlib 

matstanlib_options

%% parse required inputs
if isstruct(diagnostics)
    %actually the diagnostics structure
    if ~isfield(diagnostics,'divergent__')
        error('if input is struct-type, it must must have ''divergent__'' as a field.')
    end
    divergent = diagnostics.divergent__;
elseif isnumeric(diagnostics) && ismatrix(diagnostics) && ...
        all(ismember(diagnostics(:),[0 1]))
    %actually the divergent__ field of the diagnostics structure
    divergent = diagnostics;
else
    error(['first input must be a structure of diagnostic quantities ' ...
        '(consisitent with the output of extractsamples.m) ' ...
        'or a matrix of divergent transition indicators.'])
end

%% rug plot
nIterations = size(divergent,1);
nChains = size(divergent,2);

%start a figure ...
dumf = figure(999); %dummy figure to protect sizing
f = figure('color',[1 1 1]);
fpos = f.Position;
f.Position = [fpos([1 2]) [600 50+20*nChains]*figScaling];
close(dumf.Number); %close dummy figure
%... and an axis
axes('position',[0.125 0.1 0.775 0.3+min(nChains/10,0.45)])
hold on
axis ij %plot top to bottom

%rugplot for each chain
leftLabels = cell([nChains 1]);
rightLabels = cell([nChains 1]);
nDivByChain = sum(divergent,1);
for m = 1:nChains
    if nDivByChain(m) > 0
        divIters = find(divergent(:,m)==1);
        try
            %vertical line marker is available in R2020b or later
            plot(divIters,m*ones(size(divIters)),'|','markersize',markSz, ...
                'linewidth',linePt,'color',getcolors('red'))
        catch
            plot(divIters,m*ones(size(divIters)),'x','markersize',markSz, ...
                'linewidth',linePt,'color',getcolors('red'))
        end
    end
    leftLabels{m} = sprintf('chain %i',m);
    rightLabels{m} = num2str(nDivByChain(m));
end
%format
title('divergent transitions')
xlabel('iteration')
xlim([0 nIterations + 1])
ylim([0 nChains + 1])
%left y axis
set(gca,'ylim',[0 nChains + 1],'ytick',1:nChains,'yticklabels',leftLabels)
set(gca,'ygrid','on')
%right y axis
yyaxis right, axis ij
set(gca,'ylim',[0 nChains + 1],'ytick',1:nChains,'yticklabels',rightLabels)
set(gca,'xcolor','w','ycolor','k')
set(gca,'fontsize',fontSz)
box on

end