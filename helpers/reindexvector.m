function [outputVec,mapping] = reindexvector(inputVec,oldIndices,newIndices)
% REINDEXVECTOR reindexes a vector. 
% 
% [OUTPUTVEC,MAPPING] = REINDEXVECTOR(INPUTVEC)
%   this function will takes each unqiue element in INPUTVEC, and replaces
%   it with a sequential index.  (these elements are assigned a number
%   corresponding to their **order of appearance** in INPUTVEC.) the
%   resultant vector is returned as  OUTPUTVEC.  
%   INPUTVEC must be a vector that is either numeric or a cell of strings.
%   MAPPING, a table recounting the applied mapping from INPUTVEC elements
%   to sequential numeric indices, is also returned.  
% 
% [OUTPUTVEC,MAPPING] = REINDEXVECTOR(INPUTVEC,OLDINDICES,NEWINDICES)
%   rather than using the default mapping from unique elements to
%   sequential numbers, one can use a custom mapping, from elements of
%   OLDINDICES to corresponding elements of NEWINDICES.
%   INPUTVEC and OLDINDICES must either be both numeric type or both cells
%   of strings.  NEWINDICES must be numeric type only.
% 
%   NOTE: if there are entries in INPUTVEC that are *not found* in
%   oldIndices, then corresponding entries in OUTPUTVEC will be NaN !
% 
% 
% (c) beth baribault 2020 ---                                 > matstanlib 

%% parse inputs
if nargin == 1
    if ~(isnumeric(inputVec) || iscell(inputVec)) || ~isvector(inputVec)
        error('inputVec must be a vector of numeric or cell type.')
    end
    oldIndices = unique(inputVec,'stable');
    newIndices = 1:length(unique(inputVec));
elseif nargin == 3
    if ~(isnumeric(inputVec) || iscell(inputVec)) || ~isvector(inputVec)
        error('inputVec must be a vector of numeric or cell type.')
    elseif ~(isnumeric(oldIndices) || iscell(oldIndices)) || ~isvector(oldIndices)
        error('oldIndices must be a vector of numeric or cell type.')
    elseif ~isnumeric(newIndices) || ~isvector(newIndices)
        error('newIndices must be a vector of numeric type ONLY.')
    end
    if ~isequal(length(oldIndices),length(oldIndices))
        error('oldIndices and newIndices must be vectors of the same length.')
    end
    if any(~ismember(oldIndices,inputVec))
        error('some oldIndices are not present in inputVec, and so cannot be relabeled.')
    end
else
    error('the reindexColumn function can only be called with 1 input OR 3 inputs.')
end


%% reindex 
if isrow(oldIndices), oldIndices = oldIndices'; end %ensure column vector
if isrow(newIndices), newIndices = newIndices'; end %ensure column vector
mapping = table(oldIndices,newIndices);
if iscell(newIndices)
    outputVec = cell(size(inputVec));
elseif isnumeric(newIndices)
    outputVec = NaN(size(inputVec));
end
for n = 1:length(oldIndices)
    outputVec(ismember(inputVec,oldIndices(n))) = newIndices(n);
end

end