function [colorRGB,colorNames] = getcolors(varargin)
%GETCOLORS returns RGB values for a named color, colors, or style.
% 
% named colors are defined in colorlibrary.m.  
% 
% named styles are defined here (as lists of named colors). styles include
% color families (e.g., 'rainbow', 'muted', 'reds') and a few sets of
% colors that i like together (e.g., 'pastel'). all valid style name
% strings are listed in the first code line of this function. 
% 
% NOTE: this function is dependent on colorlibrary.m!
% NOTE: color names override style names (but that shouldn't happen).
% 
% [COLORRGB] = GETCOLORS
% [COLORRGB,COLORNAMES] = GETCOLORS
%   if the function is run with no inputs, all valid colors in the library 
%   will be returned. specifically, COLORRGB, a matrix of RBG values for 
%   each color, and COLORNAMES, a cell of color name strings, are returned.
% 
% [COLORRGB,COLORNAMES] = GETCOLORS(STRING)
% [COLORRGB,COLORNAMES] = GETCOLORS(STRING1,STRING2,...)
%   if inputs are included, they must all be of character-type.  
%   valid inputs include:
%     - a string that is a recognized color name,
%     - multiple strings that are all recognized color names,
%     - a string that is a recognized style, 
% 
% [COLORRGB,COLORNAMES] = GETCOLORS(...,'display')
%   adding 'display' as an optional input triggers the generation of a 
%   figure depicting up to 20 requested color(s).  
% 
% [COLORRGB,COLORNAMES] = GETCOLORS(...,'blend')
%   adding 'blend' as an optional input causes all colors listed in the
%   input (or found for a given style) to be blended into a single color.
% 
% NOTE: the order of the inputs to this function does not matter. 
% 
% Examples
%   [COLORRGB,COLORNAMES] = GETCOLORS('lightblue','blue','darkblue')
%   [COLORRGB,COLORNAMES] = GETCOLORS('display','rose','peach')
%   [COLORRGB,COLORNAMES] = GETCOLORS('rainbow','display')
%   COLORRGB = GETCOLORS('red','plum','white','blend','display')
% 
% See also COLORLIBRARY, MAKECOLORMAP
% 
% (c) beth baribault 2019 ---                                 > matstanlib

styleList = {'mcmc', ...
             'rainbow','base','primary','roygbiv' ...
             'dark','light','muted', ...
             'mono','monochrome','grays', ...
             'reds','oranges','yellows','greens','blues','purples', ...
             'pastel','pastels','pastels','varpastels', ...
             'ocean','oceanfade'};
optionsList = {'display','blend'};
[rgbLibrary,nameLibrary] = colorlibrary;

%parse inputs
if nargin == 0
    %if no inputs are given, return the entire color library
    colors = nameLibrary;
    blendMode = false;
    displayMode = false;
elseif nargin == 1 && isequal(varargin,{'display'})
    %if only the 'display' option is given, throw error
    error(['to use the ''display'' input option, at least one ' ...
        'color name or a style name must be supplied. ' ...
        'the entire color library cannot be displayed.'])
elseif all(ismember(varargin,optionsList))
    %if only the 'blend' option is given, throw error
    %if only the 'blend' and 'display' options are given, throw error
    error(['at least one additional input (e.g., a color name or ' ...
        'a style name) must be supplied.'])
else
    %if one or more color/style inputs are given ...
    %... first, check type of all inputs
    if ~all(cellfun(@ischar,varargin))
        error('all inputs must be strings.')
    end
    %... then, loop over inputs
    blendMode = false;    %(by default)
    displayMode = false;  %(by default)
    for v = 1:length(varargin)
        %is it a color name?
        if ismember(varargin{v},nameLibrary)
            if ~exist('colors','var')
                colors = {varargin{v}};
            else
                colors{end+1} = varargin{v};
            end
        %is it a predefined style?
        elseif ismember(varargin{v},styleList)
            style = varargin{v};
        %is it a request to blend the color(s)?
        elseif strcmp(varargin{v},'blend')
            blendMode = true;
        %is it a request to display the color(s)/style?
        elseif strcmp(varargin{v},'display')
            displayMode = true;
        else
        %is it none of the above?
            error('input string ''%s'' is not recognized.',varargin{v})
        end
    end
    if blendMode && length(colors)==1
        error('cannot blend a single color!')
    end
end
if ~exist('style','var')
    style = '   ';
end

%set color map
switch style
    %MCMC CHAINS
    case 'mcmc' %can handle up to 20 chains
        colorNames = {'darkgreen','lightgreen','darkblue','lightblue', ...
            'darkpurple','lightpurple', 'darkred','lightred', ...
            'darkorange','lightorange','darkyellow','lightyellow', ...
            'black','gray', ...
            'green','blue','purple','red','orange','yellow'};
    %core colors
    case {'rainbow','base','primary'}
        colorNames = {'red','orange','yellow','green','blue','purple'};
    case 'roygbiv'
        colorNames = {'red','orange','yellow','green','blue','indigo','purple'};
    case 'dark'
        colorNames = {'darkred','darkorange','darkyellow', ...
            'darkgreen','darkblue','darkpurple'};
    case 'light'
        colorNames = {'lightred','lightorange','lightyellow', ...
            'lightgreen','lightblue','lightpurple'};
    case 'muted'
        colorNames = {'mutedred','mutedorange','mutedyellow', ...
            'mutedgreen','mutedblue','mutedpurple'};
    case {'mono','monochrome','grays'}
        colorNames = {'black','darkgray','gray','lightgray', ...
            'lightestgray','white'};
    %color families
    case 'reds'
        colorNames = {'lightred','red','darkred','mutedred', ...
            'pink','rose','puce','lightburntorange'};    
    case 'oranges'
        colorNames = {'lightorange','orange','darkorange','mutedorange', ...
            'peach','lightpeach','mandarin','lightmandarin', ...
            'burntorange','lightburntorange','coral','lightcoral'};
    case 'yellows'
        colorNames = {'lightyellow','yellow','darkyellow','mutedyellow', ...
            'grass','chartreuse','lightchartreuse','lightpeach'};
    case 'greens'
        colorNames = {'lightgreen','green','darkgreen','mutedgreen', ...
            'chartreuse','lightchartreuse','bluegreen', ...
            'lightbluegreen','teal','mutedteal','grass'};
    case 'blues'
        colorNames = {'lightblue','blue','darkblue','mutedblue', ...
            'skyblue','teal','mutedteal','lightslate','slate','darkslate', ...
            'lightindigo'};
    case 'purples'
        colorNames = {'lightpurple','purple','darkpurple','mutedpurple', ...
            'indigo','lightindigo','plum','lightplum','lilac','fuchsia', ...
            'heliotrope'};
    %stylistic families
    case {'pastel','pastels'}
        colorNames = {'lilac','pink','peach','lightyellow', ...
            'lightgreen','mutedblue','lightplum','rose'};
    case {'varpastel','varpastels'}
        colorNames = {'lightred','lightburntorange','lightpeach','grass', ...
            'lightbluegreen','lightindigo','puce'};
    case 'ocean'
        colorNames = {'darkblue','mutedteal','lightblue', ...
            'darkgreen','mutedblue','lightslate', ...
            'bluegreen','mutedgreen','lightchartreuse', ...
            'chartreuse','teal','lightgreen','lightbluegreen', ...
            'green'};
    case 'oceanfade'
        colorNames = {'darkblue','bluegreen','darkgreen','chartreuse', ...
            'lightslate','mutedblue','mutedteal','mutedgreen', ...
            'lightblue','lightbluegreen','lightgreen','lightchartreuse'};
    otherwise
        colorNames = colors;
end
notRecognized = ~ismember(colorNames,nameLibrary);
if any(notRecognized)
    colorNames(notRecognized)
    error(['the following color names within the ''%s'' palette are ' ...
        'not recognized: %s'],style, ...
        ['''' strjoin(colorNames(notRecognized),''', ''') ''''])
end
colorRGB = rgbLibrary( ...
    cellfun(@(x) find(ismember(nameLibrary,x)),colorNames),:);

%if requested, blend colors
if blendMode
    colorRGB = mean(colorRGB,1);
    colorNames = {'< blend >'};
end

%if requested, display colors
if displayMode 
    if size(colorRGB,1) > 20
        error('only <=20 colors may be displayed at once.')
    else
        f = figure('color',[1 1 1]);
        set(f,'defaultAxesColorOrder',[0 0 0; 0 0 0])
        ax = axes('position',[0.1 0.1 0.625 0.8150]);
        hold on;
        nc = size(colorRGB,1);
        for c = 1:nc
            rectangle('pos',[0 -0.25+c 1 0.5], ...
                'facecolor',colorRGB(c,:),'edgecolor',colorRGB(c,:))
        end
        %left y axis
        set(ax,'ylim',[0.25,nc+0.75],'ytick',1:nc,'ydir','reverse')
        %right y axis
        yyaxis right
        set(ax,'ylim',[0.25,nc+0.75],'ytick',1:nc,'ydir','reverse')
        set(ax,'yticklabels',colorNames)
%         set(ax,'xlim',[0 1],'xtick',[])
        set(gca,'XColor','none','tickdir','out')
        title(style)
    end
end

end