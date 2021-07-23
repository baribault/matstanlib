function [parameter,ind] = str2ind(instanceName)
%STR2IND splits parameter instance names.
% 
% *** NOTE: this is a helper function for matstanlib, and is unlikely ***
% *** to be used directly.                                            ***
% 
% 
% [PARAMETER,IND] = STR2IND(INSTANCENAME)
%   this function will split a string, INSTANCENAME, that is a 
%   parameter instance name string of the form 'parameter[m,n,...,z]', into
%   a parameter name string, PARAMETER, and a cell of integer indices, IND.
%   
%   the indices are returned in a cell so that they might be used to
%   index into a matrix of variable size (e.g., chains(:,:,ind{:}).  
% 
% Example
%   [parameter,ind] = str2ind('lambda[3,2]');
%   lambda_3_2_chains = samples.(parameter)(:,:,ind{:});
% 
% 
% See also EXTRACTSAMPLES, EXPANDPARAMNAMES
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%check inputs
if ~ischar(instanceName)
    error('input must be a string.')
end

%extract
delims = find(instanceName=='[' | instanceName==',' | instanceName==']');
if isempty(delims)
    parameter = instanceName;
    ind = [];
else
    parameter = instanceName(1:delims(1)-1);
    dimensions = length(delims) - 1;
    if isempty(dimensions) || isequal(dimensions,0)
        ind = {1};
    else
        ind = cell([1 dimensions]);
        for d = 1:dimensions
            ind{d} = str2double(instanceName(delims(d)+1:delims(d+1)-1));
        end
    end
end

end
