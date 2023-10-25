function PerfMetrics = CalcMetrics(x,d,y,fs)

% This function computes the performance metrics ISNR, RMSE and NCC for the
% filtered signal (y) obtained through filtering algorithms from the noisy
% input signal (x) with respect to a referedence/desired signal (d).
% INPUTS  -> x (input noisy signal)                        (Nx1 vector)
%            d (desired/reference signal)                  (Nx1 vector)
%            y (filtered output)                           (Nx1 vector)
%            fs (sampling frequency)                       (scalar)
% OUTPUTS -> PerfMetrics (RMSE,ISNR,NCC metrics)           (Struct)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================
  
%% Compute ISNR (Improvement in Signal-to-Noise Ratio)
% ISNR = 20*log10(MSE(d-x)) - 20*log10(MSE(d-y)); 
% Formula taken e.g. from the book "Digital Image Processing" by S. 
% Jayaraman et al. 2011, p.357 or from the book "Encyclopedia of Optical 
% Engineering" by R.G. Driggers 2003, p.413); "d" is the original 
% image (without noise), "x" the corrupted image with noise, and 
% "y" the de-noised image.

SNR_in = norm(d - x,2); % Input SNR
SNR_out = norm(d - y,2); % Output SNR 

PerfMetrics.ISNR = 20*log10(SNR_in/SNR_out);  


%% Compute RMSE (Root Mean Squared Error)

PerfMetrics.RMSE = sqrt(mse(d - y));

%% Complex Pearson Correlation Coefficient (PCC)
% This is a value between [-1,1], calculated for the filtered and the 
% reference signal; the value -1 indicates strong anti-correlation between
% the two signals, and a value of 1 represents strong correlation. It is
% desired that PCC is close to 1.

PerfMetrics.PCC = (d'*y)/(norm(d)*norm(y));

%% R^2 (R squared)
% R^2 is the ratio of the power of the ocular artifacts being removed from 
% the primary signal to the power in the estimated EEG (the higher 
% (bounded) the value of R^2, the better the artifact minimization) [1].

% [1] S. Puthusserypady and T. Ratnarajah, "H/sup /spl infin// adaptive 
% filters for eye blink artifact minimization from electroencephalogram," 
% in IEEE Signal Processing Letters, vol. 12, no. 12, pp. 816-819, Dec. 
% 2005, doi: 10.1109/LSP.2005.859526.

PerfMetrics.R2 = norm(d - y)^2 / norm(y)^2;


%% PSD Distortion (PSD_dis)
% It calculates the distortion introduced in the PSD of the denoised signal
% w.r.t the reerence signal. It is expressed as a percentage.
% The value must be high if the algorithm works well. [1,2]

% [1] R. Ranjan, B. C. Sahana and A. K. Bhandari, "Motion Artifacts 
% Suppression From EEG Signals Using an Adaptive Signal Denoising Method," 
% in IEEE Transactions on Instrumentation and Measurement, vol. 71, pp. 
% 1-10, 2022, Art no. 4000410, doi: 10.1109/TIM.2022.3142037.
% [2] M. K. Islam, A. Rastegarnia, and Z. Yang, “A wavelet-based arti- fact
% reduction from scalp EEG for epileptic seizure detection,” IEEE J. Biomed.
% Health Inform., vol. 20, no. 5, pp. 1321–1332, Sep. 2016.

PSD_d = periodogram(d) ;
PSD_y = periodogram(y);

PerfMetrics.PSD_dis = (norm(PSD_y(1:fs/2),2)^2 / norm(PSD_d(1:fs/2),2)^2) * 100;


%% SNR_out (Output SNR in dB)
% It calculates the output SNR in dB for a given input SNR value to
% quantify the efficacy of the algorithm, as defined in [1].

% [1] Puthusserypady, S., & Ratnarajah, T. (2005). H/sup/spl infin//adaptive
% filters for eye blink artifact minimization from electroencephalogram.
% IEEE signal processing letters, 12(12), 816-819.

SNR_ref = norm(d,2)^2;
SNR_imp = norm(d-y,2)^2;

PerfMetrics.SNRout = 10*log(SNR_ref/SNR_imp);

end
    
  