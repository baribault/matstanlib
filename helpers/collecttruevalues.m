function trueValues = collecttruevalues(varargin)
%COLLECTTRUEVALUES creates a structure of true parameter values.
% 
% NOTE: this function is lazy.
% 
% TRUEVALUES = COLLECTTRUEVALUES(...)
%   this function takes all inputs and collects them in a single structure,
%   TRUEVALUES.  the added fieldnames are the same as the names of the 
%   input variables.  the added fields' values are the same as the numeric 
%   values of the input variables. 
%  
% NOTE: for consistency with other functions, row vectors will be
% transposed to become column vectors.  
% 
% Example:
%   mu = [0.1 0.5];
%   sigma = 1;
%   tv = collecttruevalues(sigma,mu);
%   tv
%   >>  struct with fields:
%   >>      sigma: 1
%   >>         mu: [2x1 double]
% 
% TRUEVALUES = COLLECTTRUEVALUES(TRUEVALUES,...)
%   this function adds all input variables' names and values to a
%   preexisting TRUEVALUES structure in the same fashion as above. 
% 
% Example (continued):
%   delta = 0.4;
%   tv = collecttruevalues(tv,delta);
%   tv
%   >>  struct with fields:
%   >>      sigma: 1
%   >>         mu: [2x1 double]
%   >>      delta: 0.4
% 
% See also GETSAMPLESTATS, PLOTRECOVERY
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% check inputs
isNumeric = cellfun(@isnumeric,varargin);
isStructure = cellfun(@isstruct,varargin);
if isStructure(1)
    if length(varargin) > 1 && all(isNumeric(2:end))
        if all(structfun(@isnumeric,varargin{1}))
            trueValuesGiven = true;
        end
    else
        error(['all inputs must be numeric type.  (the only ' ...
            'exception is if multiple inputs are given, then the ' ...
            'first input may be a structure of numeric types.)'])
    end
elseif all(isNumeric)
    trueValuesGiven = false;
else
    error(['all inputs must be numeric type.  (the only ' ...
        'exception is if multiple inputs are given, then the ' ...
        'first input may be a structure of numeric types.)'])
end

%% collect true values
%ensure any vector inputs are *column* vectors
isRowVector = cellfun(@(x) isrow(x) && ~isstruct(x),varargin);
varargin(isRowVector) = cellfun(@transpose,varargin(isRowVector),'uni',0);

%create true values structure
for n = 1:nargin
    if n==1 && trueValuesGiven
        trueValues = varargin{n};
    else
        trueValues.(inputname(n)) = varargin{n};
    end
end

end