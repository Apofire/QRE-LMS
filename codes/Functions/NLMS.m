function [y,error,wts] = NLMS(x,d,w,mu,a,ITER)
% This function computes the Normlalised LMS (NLMS) algorithm. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            w (filter coefficients/weights)               (Wx1 vector)
%            mu (step size)                                (scalar)
%            a (regularisation parameter)                  (scalar)
%            ITER (Total iterations)                       
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% INITIALISATIONS
N     = length(x);     % Input length
W     = length(w);     % Filter length
u     = zeros(W,1);    % Convolution buffer
y     = zeros(N,1);    % Filtered output
e     = zeros(N,1);    % Error per sample
error = zeros(ITER,1); % Normed error for every iteration
wts   = zeros(W,ITER); % Matrix of weights at every iteration

%% NLMS ALGORITHM
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    for i = 1:N        % Loop over signal length
        u    = [x(i);u(1:end-1,1)]; % Define signal window for convolution
        y(i) = u'*w;                % Compute filtered output by convolution
        e(i) = d(i) - y(i);         % Compute error
        w = w + ( mu/ (a + norm(u,2)^2) )*e(i)*u;  % Update the filter weights
    end
    error(iter) = norm(e,2);        % Store normed error
end
end
