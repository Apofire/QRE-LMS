function [y,err] = EKF_EEG(X,d,params,init,ITER)
% EKF based on the Dynamic State Space (DSS) equations as given in [1]

% DSS equations for EEG signal:
% for the signal length N, n = 1,2,3,..,N
% x(n) = αx(n-1) + β( x(n-1)/(1 + x^2(n-1)) ) + γcos(1.2(n-1)) + u_n
% y(n) = x^2(n)/20 + v_n
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% [1] J. Walters-Williams and Yan Li, "Comparison of Extended and 
% Unscented Kalman Filters applied to EEG signals," IEEE/ICME International
% Conference on Complex Medical Engineering, Gold Coast, QLD, Australia, 
% 2010, pp. 45-51, doi: 10.1109/ICCME.2010.5558873.
% =========================================================================

%% INITIALISATIONS 
N   = length(d);     % Length of signal
y   = init.x0;       % Initial value of the state vector
P   = init.P0*eye(N);% Initial value of the state error co-variance matrix
Q   = init.Q*eye(N); % Process noise co-variance matrix
R   = init.R*eye(N); % Measurement noise co-variance matrix
A   = diag(N);       % Linearised state transition matrix
H   = diag(N);       % Linearised measurement matrix
res = zeros(N,1);    % Residue (or innovation) vector
err = zeros(ITER,1); % Store normed error

syms x;  % Symbolic variable for state vector 
syms n;  % Symbolic variable for iteration (time) index

% DSS equations for EKF
f = params.alpha*x + params.beta*( x/(1 + x^2) ) + params.gamma*cos(1.2*n);
h = x^2/20;

% Jacobian of the DSS equations (f and h) for EKF
Ak = params.alpha + params.beta*( (1 - x^2)/(1 + x^2)^2 ); % For process map
Hk = x/10;   % For measurement map

%% EKF ALGORITHM
for iter = 1:ITER
    for i = 1:N
        A(i,i) = subs(Ak,y(i)); % Fill in the matrix A
        H(i,i)   = subs(Hk,y(i)); % Fill in the matrix H
    end
    
    % TIME UPDATE
    for i = 1:N
        % Substitute previous value of the state and the previous iteration
        % index in the process DSS equation f
        y(i) = subs(f,{x,n},{y(i),iter-1}); % Update states
    end
    P = A*P*A' + Q; % Update P

    % MEASUREMENT UPDATE
    K = P*H'/(H*P*H' + R); % Compute Kalman gain
    % Compute measurement residue (or innovation) vector
    for i = 1:N
        res(i) = X(i) - subs(h,y(i));
    end
    y = y + K*res; % State update 
    P = (eye(N) - K*H)*P; % Update P

    err(iter) = norm(d-y); % Store normed error
end


end