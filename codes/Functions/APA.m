function [y, Error, wts] = APA(x,d,Mu,ITER)
% This function implements the Affine Projection Algorithm (APA). The
% algorithm is implemented using the in-built MATLAB function for APA.

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            Mu (step size)                                (scalar)
%            ITER (Total iterations)
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% --------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% ==========================================================================

%% DEFINE PARAMETERS 
SIZE = length(x);
apfilt = dsp.AffineProjectionFilter('Length', SIZE,'StepSize', Mu);

%% INITIALISATIONS
error_apa = zeros(ITER,SIZE);
wts       = zeros(ITER,SIZE);
Error     = zeros(1,ITER);

%% ALGORITHM
for iter = 1:ITER
    [y,error_apa(iter,:)] = apfilt(x,d);     % Implement the algorithm
    wts(iter,:) = apfilt.Coefficients;       % Extract the coefficients/weights
    Error(iter) = norm(error_apa(iter,:),2); % Compute normed error
end

end