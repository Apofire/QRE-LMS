function [y,AvgError,Metrics,J] = TestFilter(Settings,fs,CleanData,...
                                           NoisyData,NoisyDataIdx,ITER,...
                                           FilterName,NumEpochs)
% This function runs an adaptive filtering algorithm for a noisy input
% signal, given a reference/clean signal, for a given number of epochs (if 
% the data is a matrix of multiple epochs), and a given number of total
% iterations. There are a total of 6 algorithms considered: LMS, logLMS,
% NLMS,ENLMS, KLMS, RLS and RENLMS. 

% Inputs  -> Settings (Filter/algorithm settings)
%              fs (Sampling frequency of the input data [for metrics calc])
%              CleanData (Clean/referecne signal)
%              NoisyData (Noise corrupted data)
%              NoisyDataIdx (Indices of the noisy data for correspondence)
%              ITER (Total number of iterations)
%              FilterName (Name of the filtering algorithm)
%              NumEpochs (Number of epochs to average over)
% Outputs -> y (Filtered output) (matrix with each column corresponding to 
%                                 one epoch)
%            AvgError (Output error averaged over 'NumEpochs' signals)
%            Metrics  (Structure storing the performance metrics) (cell array)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% INITIALISATIONS
EEGidx = NoisyDataIdx.EEG; % Collect indices used to generate noisy dataset

% Take first 'NumEpochs' columns of NoisyData
testEpochsIdx = NoisyData(:,1:NumEpochs); 

%% SWITCH BASED ON ALGORITHM
switch FilterName
    case 'LMS'   
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1); 
        Metrics    = cell(NumEpochs,1);
        for epoch  = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = LMS(x,d,Settings.w,0,ITER);     % Compute output
            Metrics{epoch}         = CalcMetrics(x,d,y(:,epoch),fs); % Store metrics
            AvgError               = AvgError + error;               % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------
    case 'IPNLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1); 
        Metrics    = cell(NumEpochs,1); 
        for epoch  = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = IPNLMS(x,d,Settings.w,Settings.mu,...
                                        Settings.alpha,ITER);        % Compute output
            Metrics{epoch}         = CalcMetrics(x,d,y(:,epoch),fs); % Store metrics
            AvgError               = AvgError + error;               % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------
    case 'NLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = NLMS(x,d,Settings.w,Settings.mu,...  
                                          Settings.a,ITER);  % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs); % Store metrics 
            AvgError = AvgError + error;                     % Compute cumulative error
        end
        AvgError = AvgError/NumEpochs;  % Compute average error
    % -------------------------------------------------------------------------    
    case 'ENLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = ENLMS(x,d,Settings.w,Settings.mu,...
                                          Settings.a,ITER);  % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs); % Store metrics 
            AvgError = AvgError + error;                     % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------    
    case 'KLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = KLMS(x,d,Settings.lambda,Settings.P,...
                                   Settings.Q,Settings.qv,Settings.w,ITER); % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs);                % Store metrics 
            AvgError = AvgError + error;                                    % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------  
    case 'KNLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = KNLMS(x,d,Settings.S_w2,Settings.qn,...
                                   Settings.qv,Settings.w,ITER);     % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs);         % Store metrics 
            AvgError = AvgError + error;                             % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------    
    case 'RLS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = RLS(x,d,Settings.M, ITER,...
                                   Settings.lambda,Settings.Rinv); % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs);       % Store metrics 
            AvgError       = AvgError + error;                     % Compute cumulative error 
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------
    case 'RENLMS'
        y          = zeros(size(CleanData,2),NumEpochs);
        AvgError   = zeros(ITER,1);
        Metrics    = cell(NumEpochs,1);
        J          = cell(NumEpochs,1);
        for epoch = 1:NumEpochs
            d = CleanData(EEGidx(epoch),:)'; % Extract clean EEG signal
            x = testEpochsIdx(:,epoch);      % Extract noisy EEG signal
            [y(:,epoch),error,wts] = RE_NLMS(x,d,Settings.w,Settings.mu,...
                                    Settings.muNLMS,Settings.a,ITER,...
                                    ITER);                     % Compute output
            Metrics{epoch} = CalcMetrics(x,d,y(:,epoch),fs);   % Store metrics 
            AvgError       = AvgError + error;                 % Compute cumulative error 
            J{epoch}       = CostFnc(wts);                     % Compute cost function
        end
        AvgError = AvgError/NumEpochs; % Compute average error
    % -------------------------------------------------------------------------   
    otherwise
        warning('Unknown Filtering Algorithm!');
        
end