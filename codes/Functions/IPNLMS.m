function [y,error,wts] = IPNLMS(x,d,w,mu,alpha,ITER)

% This function computes the Improved Proportionate NLMS (IPNLMS) as 
% proposed in [1].
% [1] F. R. Hashim, J. J. Soraghan, L. Petropoulakis and N. G. N. Daud,
% "EMG cancellation from ECG signals using modified NLMS adaptive filters,"
% 2014 IEEE Conference on Biomedical Engineering and Sciences (IECBES), 
% Kuala Lumpur, Malaysia, 2014, pp. 735-739, 
% doi: 10.1109/IECBES.2014.7047605.

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            w (filter coefficients/weights)               (Wx1 vector)
%            mu (step size)                                (scalar)
%            alpha (scaling parameter)                     (scalar)
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

%% LMS ALGORITHM
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    for i = 1:N        % Loop over signal length
        u    = [x(i);u(1:end-1,1)];          % Signal window for convolution
        y(i) = u'*w;                         % Compute output by convolution
        e(i) = d(i) - y(i);                  % Compute error
        G    = ComputeGain(W,alpha,w);       % Compute gain
        w    = w + ( mu*e(i)*G/(u'*G*u) )*u; % Update the filter weights
    end
    error(iter) = norm(e,2);        % Store normed error
end

end

function G = ComputeGain(W,alpha,w)

% This function computes the gain for the IPNLMS algorithm as given in eqs.
% (9) and (13) in [1].
% Inputs  -> W (Filter order)                         (scalar)
%            alpha (scaling parameter)                (scalar)
%            w (filter weights/coefficients)          (scalar)
% Outputs -> G (IPNLMS Gain)                          (W x W matrix)
% -------------------------------------------------------------------------

G = zeros(W,W); % Initialise output

% Compute each element of the gain
for idx = 1:W
    sumW = 0;
    for l = idx:W
        sumW = sumW + w(l);
    end
    G(idx,idx) = 0.5*( (1-alpha)/W + (1+alpha)*(abs(w(idx))/abs(sumW) ) );
end

end
