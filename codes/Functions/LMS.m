function [y,error,wts] = LMS(x,d,w,mu,ITER)
% This function computes the Least Mean Squares (LMS) algorithm. 

% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            w (filter coefficients/weights)               (Wx1 vector)
%            mu (step size)                                (scalar)
%            ITER (Total iterations)                       
% OUTPUTS -> y (filtered output)                           (Nx1 vector)
%            error (normed error)                          (ITERx1 vector)
%            wts (matrix of weights for all iterations)    (WxITER matrix)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% PREDEFINE μ IF SPECIFIED AS 0
if (mu == 0)
    
    % Find first row of correlation matrix by autocorrelating the input
    r_xx = xcorr(x);
    r_xx = r_xx(((length(r_xx)+1)/2):length(r_xx)); % Take only one-sided values
    r_xx = r_xx'; % (length(r_xx) x 1) vector

    % Form correlation matrix
    R_xx = zeros(length(r_xx),length(r_xx));
    R_xx(:,1) = r_xx;
    for col = 1:length(r_xx)-1
        r_xx = circshift(r_xx,1);
        R_xx(:,col+1) =  r_xx;
    end

    % Find maximum eigen value of the correlation matrix
    I = real(eig(R_xx));   

    mu = 2/(max(I) + min(I)); % Set the value of μ
end

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
        u    = [x(i);u(1:end-1,1)]; % Define signal window for convolution
        y(i) = u'*w;                % Compute filtered output by convolution
        e(i) = d(i) - y(i);         % Compute error
        w    = w + ( mu*e(i) )*u;   % Update the filter weights
    end
    error(iter) = norm(e,2);        % Store normed error
end

end
