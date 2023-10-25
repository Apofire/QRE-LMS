function Settings = SetSettings(M,FilterName)

% This function sets the initial parameters needed before using the
% adaptive filtering algorithms. 
% Inputs  -> M (Filter order)
%            Filtername (Name of the adaptive filtering algorithm)
% Outputs -> Settings (structure containing initial settings)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

% Check if 'Filtername' is a string
if ~isstring(FilterName)
    error('Filtername must be a string');
end

%% SWITCH BASED ON ALGORITHM
switch FilterName
    case 'LMS'
        Settings.w  = zeros(M,1);     % Filter weights
        Settings.mu = 0.0009;         % Step size
    case 'IPNLMS'
        Settings.w  = 0.1*ones(M,1);  % Filter weights
        Settings.mu = 0.001;          % Step size
        Settings.alpha = 0.5;          % Scaling parameter 
    case 'NLMS'
        Settings.w = zeros(M,1);      % Filter weights
        Settings.mu = 1.1;            % Step Size
        Settings.a  = 3;              % Regularisation parameter
    case 'ENLMS'
        Settings.w = zeros(M,1);      % Filter weights
        Settings.mu = 0.001;           % Step Size
        Settings.a  = 0.003;           % Regularisation parameter
    case 'KLMS'
        Settings.lambda = 1;          % Scaling factor for transition matrix
        Settings.P      = eye(M);     % Error co-variance matrix
        Settings.Q      = 10*eye(M);  % Process noise co-variance matrix
        Settings.qv     = 62.5;       % Measurement noise co-variance 
        Settings.w      = zeros(M,1); % Filter weights
    case 'KNLMS'
        Settings.S_w2   = 10;         % State error co-variance
        Settings.qn     = 0.1;        % Process noise  
        Settings.qv     = 0.625;      % Measurement noise 
        Settings.w      = zeros(M,1); % Filter weights 
    case 'IKLMS'
        Settings. S0_w2 = 10;         % State error co-variance
        Settings.qv     = 0.625;      % Measurement noise co-variance
        Settings.w      = zeros(M,1); % Filter weights 
    case 'RLS'
        Settings.lambda = 0.99999;    % Forgetting factor
        Settings.Rinv   = 20*eye(M);  % Inverse of auto-correlation matrix
        Settings.M      = M;          % Filter order
    case 'APA'
        Settings.mu     = 1.2;        % Step size
    case 'QRE_LMS'
        Settings.w = 2e-5*ones(M,1);% Filter weights
        Settings.mu = 0.001;         % Step Size
        Settings.a  = 3;            % Regularisation parameter
        Settings.muNLMS = 1.1;        % Step Size of NLMS (for switching)
    case 'Hinf'
        Settings.epsilon = 1.5;       % Scaling factor
        Settings.P   = 0.005*eye(M);  % Co-variance matrix
    otherwise
        warning('Algorithm does not exist!')
end







end