function varargout = interpretdiagnostics(diagnostics,rtable, ...
    maxRhat,minNeff,printReport,printWarnings)
%INTERPRETDIAGNOSTICS prints a convergence report and returns true if converged.
%   
% this function prints a report about various quantities useful for 
% assessing the quality of the mcmc sampling.  
% i.e., this function makes some human-readable statements based on 
% the stan-generated diagnostics and rhat table.)
% 
% INTERPRETDIAGNOSTICS(DIAGNOSTICS,RTABLE)
%   generates a command-line report based on a DIAGNOSTICS, a structure 
%   of directly sampled diagnostic quantities, and RTABLE, a table of 
%   calculated diagnostic quantities. 
%   by default, convergence is judged using the standard minimum
%   requirements laid out in BDA3 ...
%       -- maximum split R hat = 1.1
%       -- minimum number of effective samples = 5 per split chain
%           (i.e., 5 x the number of chains x 2)
%       -- no divergent transitions
% 
% PASSEDALLCHECKS = INTERPRETDIAGNOSTICS(...)
%   PASSEDALLCHECKS is an optional output.
%   PASSEDALLCHECKS is true if the minimum requirements were satisfied, 
%   and false otherwise.  these requirements may be  the default
%   requirements as per BDA3, or user-set critera.
% 
% INTERPRETDIAGNOSTICS(...,MAXRHAT,MINNEFF)
%   to use custom settings in determining whether PASSEDALLCHECKS is true
%   or false, use the one or both of the above listed optional inputs.
%   note that MINNEFF is the minimum *total* number of effective samples
%   for a parameter --- NOT a number to be multipled by the number of
%   chains, etc.
%   a non-zero number of divergent transitions being permissible is NOT AN
%   OPTION as divergences unequivocally indicate the corresponding
%   posterior sample is known to the results of failed sampling. 
% 
% PASSEDALLCHECKS = INTERPRETDIAGNOSTICS(DIAGNOSTICS,RTABLE,[],[], ...
%                                        PRINTREPORT,PRINTWARNINGS)
%   if PRINTREPORT is false, then no command line report will be generated.
%   if PRINTWARNINGS is false, then no warnings will be generated.
%   (the default value for both of these optional inputs is true.) 
% 
% References: Gelman, Carlin, Stern, Dunson, Vehtari, & Rubin (2013).
%               Bayesian Data Analysis, 3rd ed. CRC.
% 
% See also EXTRACTSAMPLES, RHATTABLE
% 
% (c) beth baribault 2019 ---                                 > matstanlib 

%% BDA3 defaults
defaultMaxRhat = 1.1;
defaultMinNeff = 5*2*size(diagnostics.lp__,2); %(5 per number of split chains)

%% parse inputs
if nargin < 3
    %use standard maximum Rhat (from BDA3)
    userSetMaxRhat = false;
    maxRhat = defaultMaxRhat;
else
    if isempty(maxRhat)
    %use standard maximum Rhat (from BDA3)
        userSetMaxRhat = false;
        maxRhat = defaultMaxRhat;
    else
    %use user-set maximum Rhat
        userSetMaxRhat = true;
        if ~isnumeric(maxRhat) || ~isscalar(maxRhat)
            error('maxRhat must be a single numeric value.')
        elseif maxRhat <= 1
            error('maxRhat must be greater than 1.')
        end
    end
end
if nargin < 4
    %use standard minimum number of effective samples (from BDA3)
    userSetMinNeff = false;
    minNeff = defaultMinNeff;
else
    if isempty(minNeff)
    %use standard minimum number of effective samples (from BDA3)
        userSetMinNeff = false;
        minNeff = defaultMinNeff;
    else
    %use user-set *total* minimum number of effective samples
        userSetMinNeff = true;
        if ~isnumeric(minNeff) || ~isscalar(minNeff)
            error('minNeff must be a single numeric value.')
        elseif mod(minNeff,1) > 0
            error('minNeff must be an integer.')
        elseif minNeff < size(diagnostics.lp__,2)
            error('minNeff must be greater than the number of chains.')
        end
    end
end
if nargin < 5
    printReport = true;
else
    if ~ismember(printReport,[0 1])
        error('printReport must be convertible to logical value.')
    end
end
if nargin < 6
    printWarnings = true;
else
    if ~ismember(printWarnings,[0 1])
        error('printWarnings must be convertible to logical value.')
    end
end

%% interpret diagnostics & report
vprintf(printReport,'\ninterpreting mcmc diagnostics ... \n')

%extract dimensions
nInstances = size(rtable,1);

%assess stan diagnostics
vprintf(printReport,'    the acceptance rate was:            %.3f\n', ...
    mean(diagnostics.accept_stat__(:)))
vprintf(printReport,'    the average step size was:          %.3f\n', ...
    mean(diagnostics.stepsize__(:)))
vprintf(printReport,'    the average tree depth was:         %.3f\n\n', ...
    mean(diagnostics.treedepth__(:)))

%assess divergent transitions
if any(diagnostics.divergent__(:)==1)
    percentDivergent = 100*sum(diagnostics.divergent__(:))/ ...
        numel(diagnostics.divergent__);
    if percentDivergent < 0.01
        warnStr = ['it is not advisable to use these samples as ' ...
            'the basis for inference.'];
    elseif percentDivergent < 0.5
        warnStr = 'do not use these samples as the basis for inference.';
    else
        warnStr = 'DO NOT USE THESE SAMPLES AS THE BASIS FOR INFERENCE!';
    end
    vwarning(printWarnings, ...
        ['there were %i divergent transitions (%.2f%% of samples)! ' ...
        'please use diagnostic plots to investigate potential remedies. ' ...
        '%s'],sum(diagnostics.divergent__(:)),percentDivergent,warnStr)
    vdisp(printReport,' ')
else
    vprintf(printReport,'    no divergent transitions :)\n')
end

%(don't print a table that's too long!)
maxNreport = 25;

%assess rhat values
isBadRhat = rtable.rhat > maxRhat;
if any(isBadRhat)
    if userSetMaxRhat
        %warn using user-set criterion
        vwarning(printWarnings, ...
            ['%i out of %i parameters have Rhat values above ' ...
            'the user-set criterion (> %g):'], ...
            sum(isBadRhat),nInstances,maxRhat)
    else
        %warn using standard default from BDA3
        vwarning(printWarnings, ...
            '%i out of %i parameters have bad Rhat values (> 1.1):', ...
            sum(isBadRhat),nInstances)
    end
    if sum(isBadRhat) > maxNreport
        vprintf(printReport,['\n    printing the first %i parameter instances ' ...
            'with high Rhat ...\n'],maxNreport)
        indBadRhat = find(isBadRhat);
        vdisp(printReport,rtable(indBadRhat(1:25),:))
    else
        vdisp(printReport,rtable(isBadRhat,:))
    end
else
    vprintf(printReport,'    no Rhat statistics over %g :)\n',maxRhat)
end

%assess effective sample size
isBadNeff = rtable.neff < minNeff;
if any(isBadNeff)
    if userSetMinNeff
        %warn using user-set criterion
        vwarning(printWarnings, ...
            ['the effective sample size for %i out of %i parameters ' ...
            'is below the user-set minimum effective sample size (< %i):'], ...
            sum(isBadNeff),nInstances,minNeff)
    else
        %warn using standard default from BDA3
        vwarning(printWarnings, ...
            ['the effective sample size for %i out of %i parameters ' ...
            'is too low (< 10*number of chains):'], ...
            sum(isBadNeff),nInstances)
    end
    if sum(isBadNeff) > maxNreport
        vprintf(printReport,['\n    printing the first %i parameter instances ' ...
            'with low ESS ...\n'],maxNreport)
        indBadNeff = find(isBadNeff);
        vdisp(printReport,rtable(indBadNeff(1:maxNreport),:))
    else
        vdisp(printReport,rtable(isBadNeff,:))
    end
else
    vprintf(printReport, ...
        '    no effective sample sizes under %i :)\n',minNeff)
end

vprintf(printReport,'\n')

%% determine if BDA3 (or user-set) criteria for convergence were satisfied
if any(isBadRhat) || any(isBadNeff) || any(diagnostics.divergent__(:))
    passedAllChecks = false;
else
    passedAllChecks = true;
end
if nargout == 1
    varargout{1} = passedAllChecks;
end

end

%% -------------------------------------------------------------------- %%
function vprintf(printReport,str,varargin)
    if nargin == 2 && printReport
        fprintf(str)
    elseif nargin > 2 && printReport
        fprintf(str,varargin{:})
    else
        %do not print!
    end
end

%% -------------------------------------------------------------------- %%
function vdisp(printReport,variable)
    if printReport
        disp(variable)
    else
        %do not print!
    end
end

%% -------------------------------------------------------------------- %%
function vwarning(printWarnings,str,varargin)
    if nargin == 2 && printWarnings
        warning(str)
    elseif nargin > 2 && printWarnings
        warning(str,varargin{:})
    else
        %do not print!
    end
end