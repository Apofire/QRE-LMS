function [y,error,wts] = IKLMS(x,d,S0_w2,qv,w,ITER)
% This function computes the Information Based Kalman LMS (IKLMS) 
% algorithm as propsoed in [1] (Table IV). 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            S0_w2 (State error co-variance)               (scalar)
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
PT    = 0;                 % (Eq. 32) (Truncated Co-variance)
y     = zeros(N,1);        % Filtered output 
error = zeros(ITER,1);     % Store the normed error    

%% KLMS ALGORITHM
for iter = 1:ITER % Iterate over all iterations 
    wts(:,iter) = w; % Store the weights
    for i = 1:N % Iterate over signal
        u = [x(i);u(1:end-1,1)];                     % Define signal window for convolution
        y(i) = u'*w;                                 % Filtered output
        P    = u'*u;                                 % (Eq. 20)
        alpha = d(i) - u'*w;                              % Error term 
        w = w + (conj(u)*alpha)/ (P + (PT/W) + qv/S0_w2); % Weight update
        PT = PT + P ;                                     % Co-variance update (Truncated)
    end
    error(iter) = norm(alpha,2);  % Store normed error
end

end