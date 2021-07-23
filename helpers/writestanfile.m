function varargout = writestanfile(modelCode,modelName,outputDir)
% WRITESTANFILE writes a Stan model specifiation it to a .stan file.
% 
% WRITESTANFILE(MODELCODE,MODELNAME,OUTPUTDIR)
%   this function takes the Stan model specification in MODELCODE, and
%   writes it to a new file named '<MODELNAME>.stan' located in OUTPUTDIR.
%   
%   MODELCODE must be a cell of strings, where each element is one line of
%   Stan model specification code. 
% 
% WRITESTANFILE(MODELCODE,MODELNAME)
%   if no OUTPUTDIR is given, the .stan file will be created in the current
%   working directory.  
% 
% WRITESTANFILE(MODELCODE)
%   if no MODELNAME is given, a default file name of 'untitled_model.stan'
%   will be used.
% 
% 
% STANFILEPATH = WRITESTANFILE(...)
%   optionally, the full, absolute path to the created .stan file
%   (including the file name!) may be returned as STANFILEPATH.
% 
% 
% (c) beth baribault 2021 ---                                 > matstanlib

%% check inputs
%modelCode
if ~iscell(modelCode) || ~isvector(modelCode) || ...
        ~all(cellfun(@ischar,modelCode))
    error(['the first input must be modelCode, a Nx1 cell of strings ' ...
        'containing a Stan model specification.'])
end

%%% optional inputs %%%
%modelName
if nargin < 2 || isempty(modelName)
    modelName = 'untitled_model'; %default
elseif ~ischar(modelName)
    error('modelName must be a string.')
end
%outputDir
if nargin < 3 || isempty(outputDir)
    outputDir = pwd; %default
elseif ~ischar(outputDir)
    error('outputDir must be a string representing a valid target directory.')
end

%%% outputs %%%
if nargout > 1
    error('too many outputs.')
end

%% write to .stan file
%ensure outputDir exists
if ~isfolder(outputDir)
    try
        mkdir(outputDir) %if it does not already exist, create it
    catch
        error(['failed to create the directory specified by ' ...
            'the outputDir input (%s).'],outputDir)
    end
end
%ensure correct file separator is being used for this operating system
if (isequal(filesep,'/') && contains(outputDir,'\')) || ...
        (isequal(filesep,'\') && contains(outputDir,'/'))
    badFilesep = setdiff({'/','\'},filesep); badFilesep = badFilesep{1};
    error(['outputDir contains ''%s'', but the current operating system ' ...
        'uses ''%s'' as a file separator.'],badFilesep,filesep)
end
%ensure outputDir ends in a file separator
if ~isequal(outputDir(end),filesep), outputDir(end+1) = filesep; end

%determine .stan file name & path
if isequal(modelName(end-5+1:end),'.stan')
    stanFileName = modelName; %already has .stan extension
else
    stanFileName = [modelName '.stan'];
end
stanFilePath = [outputDir stanFileName];
if isfile(stanFilePath)
    warning('overwriting current ''%s'' file in %s',stanFileName,outputDir)
end

%write the model code to a .stan file
stanFileID = fopen(stanFilePath,'w');
for n = 1:length(modelCode)
    if n < length(modelCode)
        fprintf(stanFileID,'%s\n',modelCode{n});
    else
        fprintf(stanFileID,'%s',modelCode{n});
    end
end
fclose(stanFileID);

%% return absolute path to file
if nargout==1
    varargout{1} = stanFilePath;
end

end
