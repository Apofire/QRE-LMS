function [Dataset,Indices] = Gen_DataSet(cleanSig,noiseSig,NumEpochs,SNRvals)

% This function prepares a set of data that contaminates a clean  
% signals with noise samples to achieve a specified input SNR. This is done
% using the equation x = d + λn, where x is the corrupted/contaminated
% signal, d is the clean signal, n is the noise signal, and λ is
% a scaling factor to adjust the contamination to achieve the desired input
% SNR. 
% Inputs  -> cleanSig (matrix of clean signals)
%            noiseSig (matrix of noise signals)
%            NumEpochs (number of epochs/signals to generate)
%            SNRvals (vector of different desired input SNR values)
% Outputs -> Dataset (Cell array storing contaminated signals for specific 
%                     input SNR.
%                     Each cell corresponds to a specific input SNR, 
%                     containing 'Numepochs' number of contaminated signals, 
%                     that are arranged column wise, i.e., 
%                     Dataset{SNR_val}(:,NumEpochs) is the arrangement )
%            Indices (Structure storing the randomly generated indices from 
%                     the matrix of clean and noisy signals. It contains 
%                     'NumEpochs' number of indices)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% EXTRACT DIMENSIONS AND GENERATE INDICES TO PREPARE DATASET
[EEGrows,~] = size(cleanSig); % Dimensions of clean signal matrix
[EOGrows,~] = size(noiseSig); % Dimensions of noise signal matrix

randCleanIdx = randi([1,EEGrows],[NumEpochs,1]); % Generate vector of random indices for clean signal
randNoiseIdx = randi([1,EOGrows],[NumEpochs,1]); % Generate vector of random indices for noise signal

Dataset = cell(length(SNRvals),1); % Preallocate output cell array

%% FORM THE DATASET
for cellidx = 1:length(SNRvals)
    for sigidx = 1:NumEpochs
        d = cleanSig(randCleanIdx(sigidx),:)'; % Clean EEG signal
        n = noiseSig(randNoiseIdx(sigidx),:)'; % Noise (EOG signal)
        SNRval = SNRvals(cellidx);             % Desired Input SNR
        lambda = 10^( log10( mse(d)/mse(n) ) - (SNRval)/10 ); % Compute λ
        Dataset{cellidx}(:,sigidx) = d + sqrt(lambda)*n;      % Form contaminated signal 
    end
end

% Store the indices to compare the results with the correct clean signal
Indices.EEG = randCleanIdx;
Indices.EOG = randNoiseIdx;

end
