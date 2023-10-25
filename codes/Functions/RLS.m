function [y,error,wts] = RLS(x,d,M,ITER,lambda,Rinv)
% This function computes the Recursive Least Squares (RLS) algorithm. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            M (filter order)                              (scalar)
%            ITER (Total iterations)     
%            lambda (forgetting factor)                    (scalar)
%            Rinv (inverse of auto-correlation matrix)     (MxM matrix)                 
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

% If lambda and Rinv are not specified, initialise them here
if nargin < 5 
    lambda = 0.999;
    Rinv   = 100*eye(M);
end

%% INITIALISATIONS
w     = zeros(M,1);    % Initialise filter weights
N     = length(x);     % Length of input signal
u     = zeros(M,1);    % Convolution buffer
y     = zeros(N,1);    % Filtered output
e     = zeros(N,1);    % Error per sample
error = zeros(ITER,1); % Normed error for every iteration
wts   = zeros(M,ITER); % Matrix of weights at every iteration

%% RLS ALGORITHM
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    for i = 1:N        % Loop over signal length
        u    = [x(i);u(1:end-1,1)];                     % Define signal window for convolution
        y(i) = u'*w;                                    % Compute filtered output by convolution
        e(i) = d(i) - y(i);                             % Compute error
        K    = (Rinv*u) / (lambda + u'*Rinv*u);         % Computre Gain
        Rinv = (Rinv - inv(lambda)*K*u'*Rinv) / lambda; % Update Auto-correlation
        w = w + K*(e(i));                               % Update the filter weights
    end
    error(iter) = norm(e,2);        % Store normed error
end

end