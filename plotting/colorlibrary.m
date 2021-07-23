function [rgb,names] = colorlibrary
%COLORLIBRARY has some nice named colors that are useful for plotting.  
% 
% [RGB,NAMES] = COLORLIBRARY
%   is the only valid syntax.  this returns the entire library of colors.
% 
% RGB   is a Nx3 matrix of values. each row is a RGB-01 color triplet.  
% NAMES is a NX1 cell.  each entry is a color name as a string.  
% the nth entry in NAMES corresponds to the nth row in RGB. 
% 
% See also GETCOLORS
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%base colors
primarymap(1,:)  = [.80 .12 .15]; primarynames{1}  = 'red';
primarymap(2,:)  = [.90 .44 .15]; primarynames{2}  = 'orange';
primarymap(3,:)  = [.95 .75 .22]; primarynames{3}  = 'yellow';
primarymap(4,:)  = [.45 .67 .18]; primarynames{4}  = 'green';
primarymap(5,:)  = [.00 .45 .75]; primarynames{5}  = 'blue';
primarymap(6,:)  = [.60 .25 .66]; primarynames{6}  = 'purple';
nprimaries = size(primarymap,1);

%MASTER LIST OF COLOR NAMES AND RGB-01 VALUES
%primary colors
rgb = primarymap;
names = primarynames';
%grayscale
rgb(end+1,:) = [0 0 0];      names{end+1} = 'black';
rgb(end+1,:) = 0.25*[1 1 1]; names{end+1} = 'darkgray';
rgb(end+1,:) = 0.50*[1 1 1]; names{end+1} = 'gray';
rgb(end+1,:) = 0.75*[1 1 1]; names{end+1} = 'lightgray';
rgb(end+1,:) = 0.9*[1 1 1]; names{end+1} = 'lightestgray';
rgb(end+1,:) = [1 1 1];      names{end+1} = 'white';
%primary colors, but darker
rgb(end+1:end+nprimaries,:) = primarymap*0.7;
names(end+1:end+nprimaries) = ...
    cellfun(@(x) ['dark' x],primarynames,'uni',0);
%primary colors, but lighter
rgb(end+1:end+nprimaries,:) = primarymap(1:6,:)*0.4 + 0.6;
names(end+1:end+nprimaries) = ...
    cellfun(@(x) ['light' x],primarynames,'uni',0);
%primary colors, but muted
rgb(end+1:end+nprimaries,:) = primarymap(1:6,:).*0.35 + 0.5;
names(end+1:end+nprimaries) = ...
    cellfun(@(x) ['muted' x],primarynames,'uni',0);

%other custom colors
rgb(end+1,:) = blendcolors('red','orange'); names{end+1} = 'burntorange';
rgb(end+1,:) = blendcolors('orange','yellow'); names{end+1} = 'mandarin';
rgb(end+1,:) = blendcolors('yellow','green'); names{end+1} = 'chartreuse';
rgb(end+1,:) = blendcolors('green','blue'); names{end+1} = 'bluegreen';
rgb(end+1,:) = blendcolors('blue','purple'); names{end+1} = 'indigo';
rgb(end+1,:) = blendcolors('purple','red'); names{end+1} = 'darkfuchsia';

rgb(end+1,:) = [1. .7 .6]; names{end+1}  = 'lightcoral';
rgb(end+1,:) = blendcolors('burntorange','lightcoral'); names{end+1} = 'coral';
rgb(end+1,:) = blendcolors('mandarin','lightorange'); names{end+1} = 'peach';
rgb(end+1,:) = blendcolors('mandarin','white'); names{end+1} = 'lightmandarin';

rgb(end+1,:) = [.8 .5 .5]; names{end+1}  = 'rose';
rgb(end+1,:) = [.9 .6 .6]; names{end+1}  = 'pink';
rgb(end+1,:) = [.85 .7 .8]; names{end+1}  = 'lilac';

rgb(end+1,:) = blendcolors('chartreuse','white'); names{end+1} = 'lightchartreuse';

rgb(end+1,:) = [.7 .9 1.]; names{end+1}  = 'skyblue';
rgb(end+1,:) = blendcolors('skyblue','bluegreen'); names{end+1}  = 'teal';
rgb(end+1,:) = blendcolors('teal','gray'); names{end+1}  = 'mutedteal';

rgb(end+1,:) = blendcolors('mutedblue','darkgray'); names{end+1}  = 'darkslate';
rgb(end+1,:) = blendcolors('mutedblue','gray'); names{end+1}  = 'slate';
rgb(end+1,:) = blendcolors('mutedblue','lightgray'); names{end+1}  = 'lightslate';

rgb(end+1,:) = blendcolors('bluegreen','gray'); names{end+1}  = 'graygreen';
rgb(end+1,:) = blendcolors('graygreen','white'); names{end+1}  = 'lightgraygreen';

rgb(end+1,:) = [1. .0 .7]; names{end+1}  = 'fuchsia';
rgb(end+1,:) = [.5 .1 .4]; names{end+1}  = 'plum';
rgb(end+1,:) = blendcolors('plum','white'); names{end+1}  = 'lightplum';

rgb(end+1,:) = blendcolors('lightred','lightorange'); names{end+1} = 'lightburntorange';
rgb(end+1,:) = blendcolors('lightorange','lightyellow'); names{end+1} = 'lightpeach';
rgb(end+1,:) = blendcolors('lightyellow','lightgreen'); names{end+1} = 'grass';
rgb(end+1,:) = blendcolors('lightgreen','lightblue'); names{end+1} = 'lightbluegreen';
rgb(end+1,:) = blendcolors('lightblue','lightpurple'); names{end+1} = 'lightindigo';

rgb(end+1,:) = blendcolors('mutedblue','white'); names{end+1} = 'paleblue';

%<3
rgb(end+1,:) = blendcolors('lightpurple','plum'); names{end+1} = 'heliotrope';
rgb(end+1,:) = blendcolors('lightpurple','lightred'); names{end+1} = 'puce';

%matstanlib default color
rgb(end+1,:) = blendcolors('lightbluegreen','lightblue'); names{end+1} = 'mslblue';
rgb(end+1,:) = blendcolors('lightbluegreen','lightgreen'); names{end+1} = 'mslgreen';
rgb(end+1,:) = blendcolors('lightbluegreen'); names{end+1} = 'mslbluegreen';
% rgb(end+1,:) = blendcolors('lightgreen','lightbluegreen'); names{end+1} = 'msllight';
% rgb(end+1,:) = blendcolors('skyblue','bluegreen'); names{end+1} = 'mslclear';
% rgb(end+1,:) = blendcolors('lightgreen','darkblue'); names{end+1} = 'msldark';
rgb(end+1,:) = blendcolors('mslgreen'); names{end+1} = 'msl';
rgb(end+1,:) = blendcolors('msl','white'); names{end+1} = 'msllight';
rgb(end+1,:) = blendcolors('msl','darkblue'); names{end+1} = 'msldark';

%check for duplicates
isduplicate = cellfun(@(x) sum(ismember(names,x))>1,names);
if any(isduplicate)
    warning('there are duplicate names in the color library!!!!')
    disp('indices of duplicates = ')
    disp(find(isduplicate))
    disp('names of duplicates = ')
    disp(names(isduplicate))
    error('colorlibrary will throw an error until duplicates are removed.')
end

%-------------------------------------------------------------------------%
    function [newRGB] = blendcolors(varargin)
        %input multiple color names as string, and blendcolors will average
        %their RGB values to create a new color.  returns an RGB vector. 
        ncolors = length(varargin);
        newRGB = [0 0 0];
        for v = 1:ncolors
            newRGB = newRGB + rgb(ismember(names,varargin{v}),:);
        end
        newRGB = newRGB/ncolors;
    end
%-------------------------------------------------------------------------%
    function colorRGB = getRGB(colorname)
        %input a color name as a STRING to get a vector out,
        %or input a CELL of color name strings to get a matrix out. 
        if iscell(colorname)
            colorRGB = rgb(cellfun(@(x) find(ismember(names,x)),colorname),:);
        elseif ischar(colorname)
            colorRGB = rgb(ismember(names,colorname),:);
        end    
    end
%-------------------------------------------------------------------------%
end