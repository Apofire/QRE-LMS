function [AvgMetrics,StdMetrics] = CalcAvgMetrics(Metrics,N)
% This function computes the average values of the performance metrics for
% all the "N" epochs used while testing the algorithms. 
% Inputs  -> Metrics (Performance metrics of each epoch) (cell array)
%            N (Number of test epochs)
% Outputs -> AvgMetrics (Average performance measures over N epochs)
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% INITIALISE
AvgMetrics.ISNR      = 0;
AvgMetrics.RMSE      = 0;
AvgMetrics.PCC       = 0;
AvgMetrics.R2        = 0;
AvgMetrics.PSD_dis   = 0;
AvgMetrics.SNRout    = 0;   

StdMetrics.ISNR    = 0;
StdMetrics.RMSE    = 0;
StdMetrics.PCC     = 0;
StdMetrics.R2      = 0;
StdMetrics.PSD_dis = 0;
StdMetrics.SNRout  = 0;

stdISNR    = zeros(N,1);
stdRMSE    = zeros(N,1);
stdPCC     = zeros(N,1);
stdR2      = zeros(N,1);
stdPSD_dis = zeros(N,1);
stdSNRout  = zeros(N,1);

%% COMPUTE AVERAGE
for idx = 1:N
    AvgMetrics.ISNR    = AvgMetrics.ISNR    + Metrics{idx}.ISNR;    
    AvgMetrics.RMSE    = AvgMetrics.RMSE    + Metrics{idx}.RMSE;
    AvgMetrics.PCC     = AvgMetrics.PCC     + Metrics{idx}.PCC;
    AvgMetrics.R2      = AvgMetrics.R2      + Metrics{idx}.R2;
    AvgMetrics.PSD_dis = AvgMetrics.PSD_dis + Metrics{idx}.PSD_dis;
    AvgMetrics.SNRout  = AvgMetrics.SNRout  + Metrics{idx}.SNRout;
end
AvgMetrics.ISNR    = round(AvgMetrics.ISNR/N,2);
AvgMetrics.RMSE    = round(AvgMetrics.RMSE/N,2);
AvgMetrics.PCC     = round(AvgMetrics.PCC/N,2);
AvgMetrics.R2      = round(AvgMetrics.R2/N,2);
AvgMetrics.PSD_dis = round(AvgMetrics.PSD_dis/N,2);
AvgMetrics.SNRout  = round(AvgMetrics.SNRout/N,2);


%% COMPUTE STANDARD DEVIATION 
for idx = 1:N
    stdISNR(idx)    = Metrics{idx}.ISNR;
    stdRMSE(idx)    = Metrics{idx}.RMSE;
    stdPCC(idx)     = Metrics{idx}.PCC;
    stdR2(idx)      = Metrics{idx}.R2;
    stdPSD_dis(idx) = Metrics{idx}.PSD_dis;
    stdSNRout(idx)  = Metrics{idx}.SNRout;
end

StdMetrics.ISNR    = round(std(stdISNR),2);
StdMetrics.RMSE    = round(std(stdRMSE),2);
StdMetrics.PCC     = round(std(stdPCC),2);
StdMetrics.R2      = round(std(stdR2),2);
StdMetrics.PSD_dis = round(std(stdPSD_dis),2);
StdMetrics.SNRout  = round(std(stdSNRout),2);

end


