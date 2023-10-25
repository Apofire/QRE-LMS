function [y,error,wts] = QRE_LMS(x,d,w,mu,muNLMS,a,ITER,switch_threshold)
% This function computes the Quadratic Relative Error LMS (QRE-LMS) 
% algorithm. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            w (filter coefficients/weights)               (Wx1 vector)
%            mu (step size for RENLMS)                     (scalar)
%            muNLMS (step size for NLMS)                   (scalar)
%            a (regularisation parameter)                  (scalar)
%            ITER (Total iterations)                       
%            switch_threshold (threshold to switch from RENLMS to NLMS)
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

%% RENLMS ALGORITHM
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    if iter < switch_threshold % Run RENLMS till switching threshold
        for i = 1:N
            u    = [x(i);u(1:end-1,1)]; % Define signal window for convolution
            y(i) = u'*w;                % Compute filtered output by convolution
            e(i) = d(i) - y(i);         % Compute error
            y1(i) = (u.^2)'*(w.^2);     % Temporary variable
            %            w = w + (mu*e(i)/(a + y1(i)).* ((w*w')*x_w));
            w = w + ( mu*e(i)/(a + y1(i))* (w'*w) )*u; % Update filter weights
        end
    else  % Run NLMS once threshold is reached
        for i = 1:N        % Loop over signal length
            u    = [x(i);u(1:end-1,1)]; % Define signal window for convolution
            y(i) = u'*w;                % Compute filtered output by convolution
            e(i) = d(i) - y(i);         % Compute error
            w = w + ( muNLMS/ (a + norm(u,2)^2) )*e(i)*u;  % Update the filter weights
        end
    end
    error(iter) = norm(e,2);        % Store normed error
end

end
