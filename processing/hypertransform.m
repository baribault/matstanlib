function samples = hypertransform(samples, ...
    hyperparameters,distribution,outquantity,outparameter,functions)
%HYPERTRANSFORM takes hyperparameter samples and applies known transforms.
% 
% SAMPLES = HYPERTRANSFORM(SAMPLES,HYPERPARAMETERS, ...
%               DISTRIBUTION,OUTQUANTITY,OUTPARAMETER,FUNCTIONS)
% 
% DISTRIBUTION is a string that determines which distribution 
% 
% 
% Example:
%   if in the Stan code you had:
%       alpha ~ beta(1 + a_alpha, 1 + b_alpha);       //individual-level
%       a_alpha ~ gamma(1,1);  b_alpha ~ gamma(1,1);  //group-level
%   
%   then you could obtain the group-level mean and standard deviation with:
%       samples = hypertransform(samples,{'a_alpha','b_alpha'},'beta', ...
%           {'mean','std'},{'mu_alpha','sigma_alpha'},{@(x) x+1,@(x) x+1});
% 
% (c) beth baribault 2021 ---                                 > matstanlib 

%% check inputs

%%%% check samples
if ~isstruct(samples)
    error('the first input must be the samples structure.')
end

%%%% check hyperparameters
if any(~isfield(samples,hyperparameters))
    error(['at least one hyperparameter name was not a valid field ' ...
        'in the given samples structure.'])
end

%%%% check distribution
validDistributions = {'gamma','beta','exponential','uniform'};
if ~ischar(distribution)
    error('distribution must be a string.')
elseif ~ismember(distribution,validDistributions)
    error(['distribution input not recognized.  ' ...
        'valid distribution inputs include: %s.'], ...
        ['''' strjoin(validDistributions,''', ''') ''''])
end

%%%% check quantity
if ischar(outquantity)
    outquantity = {outquantity};
end
validQuantities = {'mean','std'};
if any(~ismember(outquantity,validQuantities))
    error(['at least one quantity was not recognized.  ' ...
        'valid quantity inputs include: %s.'], ...
        ['''' strjoin(validDistributions,''', ''') ''''])
end

%%%% check outparameters
if ischar(outparameter)
    outparameter = {outparameter};
end
if any(isfield(samples,outparameter))
    error(['at least one outparameter is already a field ' ...
        'in the given samples structure.'])
end

%%%% check fhandles
if nargin < 6 || isempty(functions)
    %default is identity functions
    functions = cell(size(hyperparameters));
    for f = 1:length(functions)
        functions{f} = @(x) x;
    end
end

%% check distributional requirements
ref.gamma.nHyper = 2;
ref.gamma.names = {'shape','rate'};
ref.gamma.min = {0, 0};
ref.gamma.max = {[],[]};

ref.beta.nHyper= 2;
ref.beta.names = {'a','b'};
ref.beta.min = {0, 0};
ref.beta.max = {[],[]};

ref.exponential.nHyper= 1;
ref.exponential.names = {'rate'};
ref.exponential.min = {0};
ref.exponential.max = {[]};

ref.uniform.nHyper= 2;
ref.uniform.names = {'a','b'};
ref.uniform.min = {[], []};
ref.uniform.max = {[],[]};

%check # of hyperparameters for this distribution
nHyper = length(hyperparameters);
if ~(nHyper == ref.(distribution).nHyper)
    error(['%s distribution requires %i hyperparameter(s), ' ...
        '(%s) but %i was/were given.'], ...
        distribution,ref.(distribution).nHyper, ...
        strjoin(ref.(distribution).names,' and '),nHyper)
end
if length(functions) ~= ref.(distribution).nHyper
    error(['%s distribution has %i hyperparameter(s), ' ...
        'but %i function handle(s) was/were given.'], ...
        distribution,ref.(distribution).nHyper,length(functions))
end


%check compatability with distribution given
for h = 1:nHyper
    %minimum value
    if ~isempty(ref.(distribution).min{h})
        if any(samples.(hyperparameters{h}) < ref.(distribution).min{h})
            error(['some values of %s are below the minimum value %g ' ...
                'for the %s hyperparameter of a %s distribution.'], ...
                hyperparameters{h},ref.(distribution).min{h}, ...
                ref.(distribution).names{h},distribution)
        end
    end
    %maximum value
    if ~isempty(ref.(distribution).max{h})
        any(samples.(hyperparameters{h}) > ref.(distribution).max{h})
        error(['some values of %s are below the minimum value %g ' ...
            'for the %s hyperparameter of a %s distribution.'], ...
                hyperparameters{h},ref.(distribution).max{h}, ...
                ref.(distribution).names{h},distribution)
    end
end

%% transform
for n = 1:length(outquantity)
    switch distribution
      %gamma distribution, characterized by shape (alpha) and rate (beta)
      case 'gamma'
        f1 = functions{1};
        f2 = functions{2};
        alpha = f1(samples.(hyperparameters{1}));
        beta = f2(samples.(hyperparameters{2}));
        switch outquantity{n}
          case 'mean'
            out = alpha./beta;
          case 'median'
            error('no simple closed form for a gamma distribution median.')
          case 'mode'
            if all(alpha(:) >= 1)
                out = (alpha - 1)./beta;
            else
                error(['the mode of a gamma distribution is only defined ' ...
                    'for shape parameter (alpha) values >= 1.'])
            end
          case 'std'
            out = sqrt(alpha)./beta;
          otherwise
            error(['the quanitity ''%s'' for a %s distribution ' ...
                'is not currently defined within this function.'], ...
                outquantity{n},ditribution)
        end

      %beta distribution, characterized by shape (a) and shape (b)
      case 'beta'
        f1 = functions{1};
        f2 = functions{2};
        A = f1(samples.(hyperparameters{1}));
        B = f2(samples.(hyperparameters{2}));
        switch outquantity{n}
          case 'mean'
            out = A./(A + B);
          case 'median'
            if all(A(:) >= 1) && all(B(:) >= 1)
                out = (A - 1/3)./(A + B - 2/3); %approximate
            elseif all(A(:) == 1) && all(B(:) == 1)
                error(['the mode of a gamma distribution is defined ' ...
                    'as any value in the interval (0,1) if both shape ' ...
                    'parameters (A, B) are equal to 1.'])
            else
                error(['the mode of a gamma distribution is bimodal ' ...
                    '(specifically, is both 0 *and* 1) if both shape ' ...
                    'parameters (A, B) are < 1.'])
            end
          case 'mode'
            if all(A(:) >= 1) && all(B(:) >= 1)
                out = (A - 1)./(A + B - 2);
            elseif all(A(:) == 1) && all(B(:) == 1)
                error(['the mode of a gamma distribution is defined ' ...
                    'as any value in the interval (0,1) if both shape ' ...
                    'parameters (A, B) are equal to 1.'])
            else
                error(['the mode of a gamma distribution is bimodal ' ...
                    '(specifically, is both 0 *and* 1) if both shape ' ...
                    'parameters (A, B) are < 1.'])
            end
          case 'std'
            out = sqrt( (A.*B)./(((A + B).^2) .* (A + B + 1)) );
          otherwise
            error(['the quanitity ''%s'' for a %s distribution ' ...
                'is not currently defined within this function.'], ...
                outquantity{n},ditribution)
        end
        
      %exponential distribution, characterized by lambda (rate)
      case 'exponential'
        f1 = functions{1};
        lambda = f1(samples.(hyperparameters{1}));
        switch outquantity{n}
          case 'mean'
            out = 1./lambda;
          case 'median'
            out = log(2)./lambda;
          case 'mode'
            out = 0;
          case 'std'
            out = sqrt(1./(lambda.^2));
          otherwise
            error(['the quanitity ''%s'' for a %s distribution ' ...
                'is not currently defined within this function.'], ...
                outquantity{n},ditribution)
        end
        
      %uniform distribution, characterized by the min (a) and max (b) values
      case 'uniform'
        f1 = functions{1};
        f2 = functions{2};
        a = f1(samples.(hyperparameters{1}));
        b = f2(samples.(hyperparameters{2}));
        switch outquantity{n}
          case 'mean'
            out = 0.5.*(a + b);
          case 'median'
            out = 0.5.*(a + b);
          case 'mode'
            error(['the mode of a continuous uniform distribution is ' ...
                'any value in (a,b).'])
          case 'std'
            out = sqrt(1/12).*(b - a);
          otherwise
            error(['the quanitity ''%s'' for a %s distribution ' ...
                'is not currently defined within this function.'], ...
                outquantity{n},ditribution)
        end
    end
    samples.(outparameter{n}) = out;
end

end
