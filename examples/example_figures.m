%matstanlib paper figures
% 
%(c) beth baribault 2020 ---

clear
clc
% close all

fprintf('\n\n*** matstanlib ***\n\n')

%% inputs
%save figures? where? 
figFolder = [pwd filesep 'figures'];
figFormats = {'pdf','png'};%{'pdf','png'};

%%
%figure saving
if figFolder(end)~=filesep, figFolder(end+1) = filesep; end
if ~isfolder(figFolder), mkdir(figFolder); end
for ff = 1:length(figFormats)
    if figFormats{ff}(1)~= '.', figFormats{ff} = ['.' figFormats{ff}]; end
%     figFormats{ff} = ['_' figFormats{ff}];
end
figure % (ensures MATLAB doesn't goof the sizing on the first plot)
title('dummy figure (to protect figure sizing)')

%% funnel
close all 
figure % (ensures MATLAB doesn't goof the sizing on the first plot)
title('dummy figure (to protect figure sizing)')

%load output data
load('output_funnel.mat')

%%%%%%%%%% joint density
jointdensity(samples,'x[1]','y',diagnostics)

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_jointdensity' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- jointdensity\n',f.Number)

%%%%%%%%%% plotautocorr
plotautocorr(samples,'x[1]')

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_plotautocorr' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- plotautocorr\n',f.Number)

%%%%%%%%%% plotautocorr, with 'trace'
plotautocorr(samples,'x[1]','trace')

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_plotautocorr_trace' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- plotautocorr, with ''trace'' \n',f.Number)

%%%%%%%%%% plotlp
plotlp(diagnostics,rtable)

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_plotlp' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- plotlp\n',f.Number)

%%%%%%%%%% plotdivergences
plotdivergences(diagnostics)

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_plotdivergences' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- plotdivergences\n',f.Number)

%%%%%%%%%% parcoordivergent
parcoordivergent(samples,diagnostics,parameterRequest)

for ff = 1:length(figFormats)
    exportgraphics(gcf,[figFolder 'funnel_parcoordivergent' figFormats{ff}])
end
%echo figure
f = gcf; fprintf('FIGURE %i: funnel --- parcoordivergent\n',f.Number)


%% RL --- bad priors
close all 
figure % (ensures MATLAB doesn't goof the sizing on the first plot)
title('dummy figure (to protect figure sizing)')

%load output data
load('output_RL_badpriors.mat')