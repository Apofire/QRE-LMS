%% MAIN SCRIPT TO IMPLEMENT ALL ALGORITHMS
% This sript implements all functions for the given datasets for EEG signal
% filtering/system identification. 

% METHODS:
% 1. Least Mean Squares (LMS)
% 2. Normalised Least Mean Squares (NLMS)
% 3. Error Normalised Least Mean Squares (ENLMS)
% 4. Recursive Least Squares (RLS) 
% 5. Kalman LMS (KLMS)
% 6. Normalised Kalman LMS (KNLMS)
% 7. Information based Kalman LMS (IKLMS)
% 8. H-infinity method (Hinf)
% 9. Improved Proportionate NLMS (IPNLMS)
% 10. Proposed QRE-LMS 
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% ADD THE DIRECTORY CONTAINING THE DATASET 

addpath('../codes/Dataset/'); % Add the folder containing the datasets
load EEGDATASET.mat; % Load EEG dataset
load EMGDATASET.mat; % Load EMG dataset
load EOGDATASET.mat; % Load EOG dataset
addpath('../codes/Functions/'); % Add folder containing functions

%% GENERATE DATASETS OR LOAD IF ALREADY GENERATED
NumEpochs = 100; % Number of signals/epochs for dataset
SNRvals   = [-2, -5, -10, -20, -25,]; % Different input SNR values in dB

% GENERATE EEG + EOG SIGNAL DATASET IF NOT GENERATED ALREADY
if ~exist('../codes/Dataset/EEG_EOG_DATASET.mat','file')
    [EEG_EOG_Dataset,EEG_EOG_indices] = Gen_DataSet(EEG_all_epochs,...
                                        EOG_all_epochs,NumEpochs,SNRvals);
    save('../codes/Dataset/EEG_EOG_DATASET.mat', 'EEG_EOG_Dataset', ...
         'EEG_EOG_indices');
else 
    load EEG_EOG_DATASET.mat;
end


% GENERATE EEG + EMG SIGNAL DATASET IF NOT GENERATED ALREADY
if ~exist('../codes/Dataset/EEG_EMG_DATASET.mat','file')
    [EEG_EMG_Dataset,EEG_EMG_indices] = Gen_DataSet(EEG_all_epochs,...
                                        EMG_all_epochs,NumEpochs,SNRvals);
    save('../codes/Dataset/EEG_EMG_DATASET.mat', 'EEG_EMG_Dataset', ...
         'EEG_EMG_indices');
else
    load EEG_EMG_DATASET.mat;
end

%% PLOT THE WELCH SPECTRUM OF A RANDOM CLEAN AND NOISY EEG SIGNAL
% 
% % Plot for EOG contamination
% EOGidx = PlotSampleSignals(EEG_all_epochs,EEG_EOG_Dataset,EEG_EOG_indices,fs,SNRvals,"Freq");
% 
% % Plot for EMG contamination
% EMGidx = PlotSampleSignals(EEG_all_epochs,EEG_EMG_Dataset,EEG_EMG_indices,fs,SNRvals,"Freq");
% 
% %% PLOT THE TIME DOMAIN PLOTS OF A RANDOM CLEAN AND NOISY EEG SIGNAL
% 
% % Plot for EOG contamination
% PlotSampleSignals(EEG_all_epochs,EEG_EOG_Dataset,EEG_EOG_indices,fs,SNRvals,"Time",EOGidx);
% 
% % Plot for EMG contamination
% PlotSampleSignals(EEG_all_epochs,EEG_EMG_Dataset,EEG_EMG_indices,fs,SNRvals,"Time",EMGidx);

%% TEST ALGORITHMS 

% INITIAL SETTINGS FOR ALGORITHMS
M = 512; % Filter order 
LMS_settings     = SetSettings(M,"LMS"); 
IPNLMS_settings  = SetSettings(M,"IPNLMS");
NLMS_settings    = SetSettings(M,"NLMS");
ENLMS_settings   = SetSettings(M,"ENLMS");
KLMS_settings    = SetSettings(M,"KLMS");
KNLMS_settings   = SetSettings(M,"KNLMS");
IKLMS_settings   = SetSettings(M,"IKLMS");
RLS_settings     = SetSettings(M,"RLS");
QRE_LMS_settings = SetSettings(M,"QRE_LMS");

ITER = 1000;         % Total number of iterations
NumTestEpochs = 1; % Number of epochs for testing

% RUN THE ALGORITHMS AND STORE OUTPUTS
% [yLMS,AvgErLMS,LMSMetrics] = TestFilter(LMS_settings,fs, ...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"LMS",NumTestEpochs);
% 
% [yIPNLMS,AvgErlogIPNLMS,IPNLMSMetrics] = TestFilter(IPNLMS_settings,fs, ...
%                             EEG_all_epochs,EEG_EMG_Dataset{5},...
%                             EEG_EMG_indices,ITER,"IPNLMS",NumTestEpochs);
% 
% [yNLMS,AvgErNLMS,NLMSMetrics] = TestFilter(NLMS_settings,fs,...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"NLMS",NumTestEpochs);
% 
% [yENLMS,AvgErENLMS,ENLMSMetrics] = TestFilter(ENLMS_settings,fs,...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"ENLMS",NumTestEpochs);
% 
% [yKLMS,AvgErKLMS,KLMSMetrics] = TestFilter(KLMS_settings,fs,...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"KLMS",NumTestEpochs);
% 
% [yKNLMS,AvgErKNLMS,KNLMSMetrics] = TestFilter(KNLMS_settings,fs,...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"KNLMS",NumTestEpochs);
% 
% [yRLS,AvgErRLS,RLSMetrics] = TestFilter(RLS_settings,fs,...
%                             EEG_all_epochs,EEG_EMG_Dataset{1},...
%                             EEG_EMG_indices,ITER,"RLS",NumTestEpochs);

[yQRELMS,AvgErQRELMS,QRELMSMetrics,JQRELMS] = TestFilter(QRE_LMS_settings,fs,...
                            EEG_all_epochs,EEG_EMG_Dataset{1},...
                            EEG_EMG_indices,ITER,"QRE_LMS",NumTestEpochs);

%% COMPUTE AVERAGE METRICS

[LMSAvgMetrics, LMSStd]       = CalcAvgMetrics(LMSMetrics,NumTestEpochs);
[IPNLMSAvgMetrics, IPNLMSStd] = CalcAvgMetrics(IPNLMSMetrics,NumTestEpochs);
[NLMSAvgMetrics, NLMSStd]     = CalcAvgMetrics(NLMSMetrics,NumTestEpochs);
[ENLMSAvgMetrics, ENLMSStd]   = CalcAvgMetrics(ENLMSMetrics,NumTestEpochs);
[KLMSAvgMetrics, KLMSStd]     = CalcAvgMetrics(KLMSMetrics,NumTestEpochs);
[KNLMSAvgMetrics, KNLMSStd]   = CalcAvgMetrics(KNLMSMetrics,NumTestEpochs);
[RLSAvgMetrics, RLSStd]       = CalcAvgMetrics(RLSMetrics,NumTestEpochs);
[QRELMSAvgMetrics, QRELMSStd] = CalcAvgMetrics(QRELMSMetrics,NumTestEpochs);

%% Plot the SNRout for different values of M
% PlotSNRout();

%% COMPARE CONVERGENCE OF WEIGHTS TO OPTIMAL VALUE
OptimalWeightComparison(EEG_all_epochs,EEG_EOG_Dataset{1},EEG_EOG_indices,...
    LMS_settings,NLMS_settings,ENLMS_settings,KLMS_settings,RLS_settings,...
    QRE_LMS_settings,ITER);

