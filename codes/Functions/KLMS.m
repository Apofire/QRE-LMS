function [y,error,wts] = KLMS(x,d,lambda,P,Q,qv,w,ITER)
% This function computes the Kalman-LMS algorithm as propsoed in [1]. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            lambda (scaling factor for transition matrix) (scalar)
%            P (State error co-variance matrix)            (WxW matrix)
%            Q (Process noise co-variance matrix)          (WxW matrix)
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
L2    = lambda^2;         

%% KLMS ALGORITHM
for iter = 1:ITER % Iterate over all iterations 
    wts(:,iter) = w; % Store the weights
    for i = 1:N % Iterate over signal
        u = [x(i);u(1:end-1,1)];                            % Define signal window for convolution
        y(i) = u'*w;                                        % Filtered output
        alpha = d(i) - u'*w;                                % Error term 
        w = lambda*w + lambda*(P*u*alpha)/ (u'*P*u + qv); % Weight update
        P = L2*P - L2*((P*u)*(u'*P))/(u'*P*u + qv) + Q;     % Co-variance update
    end
    error(iter) = norm(alpha,2);  % Store normed error
end

end