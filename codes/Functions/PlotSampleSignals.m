function idx = PlotSampleSignals(EEG_all_epochs,NoisyGenDataset,DatasetIndices,Fs,SNRvals,Domain,idx)
% This function plots the frequency (Welch) spectrum or the time domain
% representation of a randomly chosen EEG signal and a noise contanimated
% EEG signal (with EMG or EOG) for a given input SNR value.
% Inputs  -> EEG_all_epochs (all clean EEG signals)
%            NoisyGenDataset (Noisy dataset generated using EOG/EMG)
%            DatasetIndices (Indices to keep track of clean EEG signals)
%            Fs (Sampling frequency)
%            SNRvals (index for "SNRvals" array to choose input SNR)
%            Domain (Type of plot: 'Freq' or 'Time')
%            idx [optional] (index of the signal to plot)
% Outputs -> Plot of Spectra (Clean and Noisy EEG)
%            idx [optional] to save the index of the signal used (in case 
%                           to compare frequency and time plots)
% -------------------------------------------------------------------------
% Code written in part by: KAUSHIK IYER 
% =========================================================================

%% EXTRACT SIGNALS
if nargin < 7 % If idx not specified as an argument
    idx       = randi(length(DatasetIndices.EEG)); % Random index to choose EEG signal
end
L         = length(SNRvals);                   % Number of input SNRs 
EEGidx    = DatasetIndices.EEG(idx);           % Choose a random index in dataset
EEGidxPos = find(DatasetIndices.EEG == EEGidx);% Find the position of EEGidx

CleanEEG  = EEG_all_epochs(EEGidx,:)';    % Clean EEG
N         = length(CleanEEG);             % Length of the signal
NoisyEEG  = zeros(N,length(SNRvals));     % Preallocate Matrix       
% Extract Noisy signal for each input SNR value
for snr = 1:L
    NoisyEEG(:,snr)  = NoisyGenDataset{snr}(:,EEGidxPos)'; % Noisy EEG
end

 %% PLOTS
if Domain == "Freq" % FREQUENCY SPECTRUM
    % Plot half of the power spectrum with 50% overlap and Chebwin window of
    % length 128
    [powerSpecClean,freqVecClean] = pwelch(CleanEEG,chebwin(128,100),[],N,Fs);
    figure;  % Define figure window
    for snr = 1:L
        [powerSpecNoisy,freqVecNoisy] = pwelch(NoisyEEG(:,snr),chebwin(128,100),[],N,Fs);
        subplot(ceil(L/2),2,snr)
        plot(freqVecClean,20*log10(abs(powerSpecClean)),'LineWidth',2,'Color','b')
        hold on
        plot(freqVecNoisy,20*log10(abs(powerSpecNoisy)),'LineWidth',2,'Color','r')
        TitleName = ['Welch spectrum for input SNR = ',num2str(SNRvals(snr)), 'dB'];
        title(TitleName,'Interpreter','latex','FontSize',24)
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',24)
        ylabel('Magnitude (dB)','Interpreter','latex','FontSize',24);
        legend('Clean Signal','Noisy Signal')

        % Axis settings
        ax = gca;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridLineStyle = '--';
        ax.FontSize = 20;
        ax.LineWidth = 2;
        box on
    end

elseif Domain == "Time" % TIME DOMAIN
    figure;  % Define figure window
    for snr = 1:L
        subplot(ceil(L/2),2,snr)
        plot(CleanEEG,'LineWidth',2,'Color','b')
        hold on
        plot(NoisyEEG(:,snr),'LineWidth',2,'Color','r')
        TitleName = ['Input SNR = ',num2str(SNRvals(snr)), 'dB'];
        title(TitleName,'Interpreter','latex','FontSize',24)
        xlabel('Samples','Interpreter','latex','FontSize',24)
        ylabel('Aplitude ($\mu V$)','Interpreter','latex','FontSize',24);
        legend('Clean Signal','Noisy Signal')

        % Axis settings
        ax = gca;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridLineStyle = '--';
        ax.FontSize = 20;
        ax.LineWidth = 2;
        box on
    end
else 
    error("Incorrect Choice of plot! Choose 'Freq' or 'Time'");
end

end