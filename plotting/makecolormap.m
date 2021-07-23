function map = makecolormap(color1,color2,nPoints)
%MAKECOLORMAP returns a colormap for a gradient between two colors.
% 
% MAP = MAKECOLORMAP(COLOR1, COLOR2)
%   both COLOR1 and COLOR2 must be 1x3 vectors representing RBG-01 
%   colors, or recognized color names (i.e., included in colorlibrary.m).  
%   this returned variable, MAP, is a 100x3 numeric colormap.
% 
% MAP = MAKECOLORMAP(...,NPOINTS)
%   NPOINTS may be used to adjust the apparent smoothness of the 
%   colormap gradient. MAP is a NPOINTSx3 numeric colormap.
% 
% Examples:
%   MAP = MAKECOLORMAP([1 0 0],[0 0 1],128)
%   MAP = MAKECOLORMAP(getcolors('red'),getcolors('blue'))
% 
% See also COLORMAP, GETCOLORS, COLORLIBRARY
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% check inputs
%color1
if ischar(color1)
    color1 = getcolors(color1);
elseif ~isnumeric(color1) || ~isequal(size(color1),[1 3]) || ~all(0 <= color1 & color1 <=1)
    error(['each color input must be a 1x3 vector ' ...
        'representing a RGB-01 color.'])
end
%color2
if ischar(color2)
    color2 = getcolors(color2);
elseif ~isnumeric(color2) || ~isequal(size(color2),[1 3]) || ~all(0 <= color2 & color2 <=1) 
    error(['each color input must be a 1x3 vector ' ...
        'representing a RGB-01 color.'])
end
if nargin < 3
    nPoints = 64; %same default as MATLAB's colormaps
end

%% colormap
map = [linspace(color1(1),color2(1),nPoints)' ...
       linspace(color1(2),color2(2),nPoints)' ...
       linspace(color1(3),color2(3),nPoints)'];

end