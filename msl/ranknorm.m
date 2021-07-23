function z = ranknorm(chains,offset)
%RANKNORM rank normalizes posterior samples with a fractional offset.  
% 
% Z = RANKNORM(CHAINS)
%   the posterior samples in the CHAINS matrix are ranked (after pooling
%   across chains), then normalized with a fractional offset of 3/8.  
% 
% Z = RANKNORM(X,OFFSET)
%   a different fractional offset may also be given.  
% 
% (c) beth baribault 2021 ---                                 > matstanlib

if nargin < 2
    %as per Vehtari et al. (2020), with reference to Blom (1958)
    offset = 3/8;
end

sz = size(chains);  %extract size
S = prod(sz);       %total number of samples
z = chains(:);      %convert from matrix to vector

%rank normalize
z = tiedrank(z);    %rank
z = (z - offset)/(S - 2*offset + 1); %fractional offset
z = norminv(z,0,1); %normalize

z = reshape(z,sz);  %convert from vector to matrix
end