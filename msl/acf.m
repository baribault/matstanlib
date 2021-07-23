function [auco,lags] = acf(x)
%ACF estimates the autocorrelation function via FFT.
% 
% [AUCO, LAGS] = ACF(X)
%   computes the biased estimate (1/N) of the autocorrelation function,
%   AUCO, for non-negative LAGS (i.e., 0:N-1), from the vector X via fast
%   fourier transform (for greater efficiency and computational accuracy). 
% 
% (c) beth baribault 2021 ---                                 > matstanlib

if ~isvector(x)
    error('only vector inputs are accepted.')
elseif ~isnumeric(x)
    error('input vector must be numeric.')
end

S = numel(x);       %total number of samples

%compute autocorrelation via FFT
%(and avoid requiring Signal Processing Toolbox just for xcorr)
nZeros = 2^nextpow2(2*S - 1);
x_fft = fft(x - mean(x),nZeros); %%%%%%%%%%%% subtract off mean
auco = ifft(x_fft.*conj(x_fft));
auco = auco/S; %(use the biased estimate, as per Vehtari et al., 2020)

%shift such that lags are in ascending order, and return as column vectors
% lags = -(S-1):(S-1); lags = lags';
% auco = [auco(end-S+2:end);auco(1:S)];
lags = (0:S-1)';
auco = auco(1:S);

end
