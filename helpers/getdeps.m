function varargout = getdeps(functionName,echoReport)
%GETDEPS echoes the function & toolbox dependencies for a function. 
% 
% MATLAB's syntax for doing this is annoying. so i wrote this. wow
% (possibly may write getalldeps.m in future to add some recursion?)
% 
% getdeps(FUNCTIONNAME)
%   will echo a report at the command line.
% 
% [FUNCTIONS,TOOLBOXES] = GETDEPS(FUNCTIONNAME)
%   will not echo a report.
% 
% NOTE: calling MATLAB's built-in dependency function from within a 
% function is desirable because it effectively clears the workspace, 
% thus avoiding ambiguity due to function handles currently in the 
% workspace, etc. 
% 
% NOTE: MATLAB's GUI-based "Dependency Reports" may also be used. 
% 
% (c) beth baribault 2019---                                  > matstanlib 

%% parse inputs
if ischar(functionName)
    if ~strcmp(functionName(end-1:end),'.m')
        functionName = [functionName '.m'];
    end
    if ~exist(functionName,'file')
        error('the function ''%s'' does not exist.',functionName)
    end
else
    error('functionName must be submitted as a string.')
end
if nargin == 1
    if nargout == 0
        echoReport = true;
    elseif nargout == 2
        echoReport = false;
    else
        error('only 0 or 2 outputs are permitted.')
    end
end

%% get dependencies
if verLessThan('matlab','8.3') %before v8.3, R2014a
    error(['this code is not backward compatible with your version' ...
        'of MATLAB.  try calling DEPFUN instead.'])
    %NOTE: depfun deprecated in v9.0, R2016a
else
    [fList,pList] = ...
        matlab.codetools.requiredFilesAndProducts(functionName);
end

%reformat for output
functions = fList';
toolboxes = {pList.Name}';

%if requested, command line report
if echoReport
    fprintf('functions required by %s:\n',functionName)
    disp(functions)
    fprintf('toolboxes required by %s:\n',functionName)
    disp(toolboxes)
    if isempty(functions) && isequal(toolboxes,{'MATLAB'})
        fprintf('... this suggests that %s is a built-in.\n',functionName)
    end
end

%% parse outputs
if nargout == 2
    varargout{1} = functions;
    varargout{2} = toolboxes;
end
    
end