function [y,error,wts] = Hinf_AF(x,d,M,ITER,epsilon,P)
% This function computes the H-infinity (Hinf) algorithm. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            M (filter order)                              (scalar)
%            ITER (Total iterations)     
%            epsilon (scaling factor)                      (scalar)
%            P (Co-variance matrix)                        (MxM matrix)                 
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

% If lambda and Rinv are not specified, initialise them here
if nargin < 5 
    epsilon = 1.5;
    P   = 0.005*eye(M);
end

%% INITIALISATIONS
w     = zeros(M,1);    % Initialise filter weights
N     = length(x);     % Length of input signal
u     = zeros(M,1);    % Convolution buffer
y     = zeros(N,1);    % Filtered output
e     = zeros(N,1);    % Error per sample
error = zeros(ITER,1); % Normed error for every iteration
wts   = zeros(M,ITER); % Matrix of weights at every iteration
eg    = (1/epsilon)^2; % Variable to be used in algorithm

%% H_infinity ALGORITHM
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    for i = 1:N        % Loop over signal length
        u      = [x(i);u(1:end-1,1)];                    % Define signal window for convolution
        y(i)   = u'*w;                                   % Compute filtered output by convolution
        e(i)   = d(i) - y(i);                            % Compute error
        r      = u*u';                                   % Correlation matrix
        P      = inv(inv(P) + (1 - eg)*r) + 1e-6*eye(M); % Co-variance update
        PTilde = inv(inv(P) - eg*r);                     % Scaled co-variance
        K      = (PTilde*u) / (1 + u'*PTilde*u);         % Computre Gain
        w      = w + K*(e(i));                           % Update the filter weights
    end
    error(iter) = norm(e,2);        % Store normed error
end

end