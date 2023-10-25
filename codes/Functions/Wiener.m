function [y,error,wts] = Wiener(x,d,W)
% This function implements the Wiener filter based on the Wiener-Hopf 
% equations.

% INPUTS  -> x (Noisy signal)
%            d (Reference/desired signal)
%            W (Filter order)
% OUTPUTS -> y (Filtered signal)
%            wts (Wiener filter weights)
%            error (normed error)
% -------------------------------------------------------------------------
% Code modified by: KAUSHIK IYER 
% =========================================================================

X = 1/W .* fft(x(1:W)); % FFT of input signal
Y = 1/W .* fft(d(1:W)); % FFT of reference signal
X = X(:);
Y = Y(:);

% Compute Power functions
Rxx = W .* real(ifft(X .* conj(X))); % Autocorrelation function
Rxy = W .* real(ifft(X .* conj(Y))); % Crosscorrelation function
Rxx = toeplitz(Rxx); % Create a circulant matrix
Rxy = Rxy';

% Wiener-Hopf vector-matrix equation to calculate weights
wts = Rxy / Rxx; wts = wts(:); 

y = fftfilt(wts,x); % FFT convolution to compute the output (filtered) signal
y = y(W+1:end); % Store only the last W terms 

error = norm((d(W+1:end) - y),2)^2; % Normed error

end