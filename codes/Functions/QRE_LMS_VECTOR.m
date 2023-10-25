%% RE-NLMS FILTER WITH VECTOR/MATRIX UPDATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Somehow takes longer in computation than the scalr counterpart in 
% RE_NLMS.m and the result is off, probably because of the diagonalisation
% of the weight matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,error,wts] = QRE_LMS_VECTOR(x,d,w,mu,muNLMS,a,ITER,switch_threshold)
% This function computes the QRE-LMS algorithm in a 
% vector/matrix format. 

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

%% RENLMS ALGORITHM (VECTOR)
for iter = 1:ITER      % Loop over iterations
    wts(:,iter) = w;   % Store the weights
    if iter < switch_threshold % Run RENLMS till switching threshold
        for i = 1:N
            u    = [x(i);u(1:end-1,1)]; % Define signal window for convolution
            D     = diag(w);                       % Form a diagonal matrix of the weights
            y(i) = u'*w;                           % Compute filtered output by convolution
            e(i) = d(i) - y(i);                    % Define the error signal
            Gamma_a = (D*u)/( norm(D*u,2)^2 + a ); % Define Gamma_a
            w = w + (mu/i)*e(i)*D'*Gamma_a;        % Update the filter weights
        end
    else        % Run NLMS once threshold is reached
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
