function [y,error,wts] = KNLMS(x,d,S_w2,qn,qv,w,ITER)
% This function computes the Normalised Kalman-LMS algorithm as propsoed 
% in [1] (Table III). 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            S_w2 (State error co-variance)                (scalar)
%            qn (Process noise co-variance)                (scalar)
%            qv (measurement noise co-variance matrix)     (scalar)
%            w (filter coefficients/weights)               (Wx1 vector)
%            ITER (Total iterations)
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% --------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% 
% [1] Lopes, P. A. C., & Gerald, J. B. (2007). New Normalized LMS Algorithms
%     Based on the Kalman Filter. 2007 IEEE International Symposium on 
%     Circuits and Systems. doi:10.1109/iscas.2007.378235
% ==========================================================================

%% Initialisations 
N     = length(x);         % Length of input vector
W     = length(w);         % Length of filter weights
u     = zeros(W,1);        % Convolution buffer
wts   = zeros(W,ITER);     % Store all the weights
alpha = zeros(N,1);        % Error per sample
y     = zeros(N,1);        % Filtered output 
error = zeros(ITER,1);     % Store the normed error    

%% KLMS ALGORITHM
for iter = 1:ITER % Iterate over all iterations 
    wts(:,iter) = w; % Store the weights
    for i = 1:N % Iterate over signal
        u = [x(i);u(1:end-1,1)];                     % Define signal window for convolution
        y(i) = u'*w;                                 % Filtered output
        P    = u'*u;                                 % (Eq. 20)
        alpha = d(i) - u'*w;                         % Error term 
        w = w + (conj(u)*alpha)/ (P + qv/S_w2);           % Weight update
        S_w2 = S_w2*(1 - (P/W)/(P + qv/S_w2)) + qn ; % Co-variance update
    end
    error(iter) = norm(alpha,2);  % Store normed error
end

end