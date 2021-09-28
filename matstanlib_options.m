%MATSTANLIB_OPTIONS
% 
% this script sets options for various matstanlib functions.  
% it is intended to be called from within matstanlib functions *only*.
% 
% (c) beth baribault 2020 ---                                 > matstanlib

%% how should figures be formatted?
%virtually all users of matstanlib should use the default option:
%       >>>  'default'
%but a very small number of users who have a hiDPI (or "pixel-doubling")
%monitor and are running an OS where MATLAB does not offer hiDPI support 
%(e.g., those running Linux on some Dell XPS laptops) should use the
%alternative option: 
%       >>>  'hiDPI'
%however, you may also elect to use custom figure formatting option values.
%PLEASE NOTE that if you engage this option, you risk breaking matstanlib.  
%       >>>  'custom'
% 
%(if you're not sure, leave this set to 'default'.) 
figureFormattingMode = 'default';

%if you have set 'custom' above and change the values below, then
%matstanlib will use your custom values:  
figScaling = 1;   %figure size (proportionally resizes the entire figure)
fontSz = 10;      %font size
markSz = 8;       %marker size
linePt = 1;       %line width






























%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%% %%%%%%  deleting this block of code will break matstanlib.  %%%%%% %%%
%%%% %%%%%%   please do not change anything below this line!   %%%%%% %%%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
  
%figure & text sizing
switch figureFormattingMode
    case 'custom'
        %use the custom settings, from above
        % ...             %figure size
        % ...             %font size
        % ...             %marker size
        % ...             %line width
    case 'hiDPI'
        %for hiDPI monitors (if figures look small & weirdly proportioned)
        figScaling = 1;   %figure size
        fontSz = 18;      %font size
        markSz = 12;      %marker size
        linePt = 1.5;     %line width
    otherwise %'default'
        %for regular monitors (the default options)
        figScaling = 1;   %figure size
        fontSz = 10;      %font size
        markSz = 8;       %marker size
        linePt = 1;       %line width
end
