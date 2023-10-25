function PlotSNRout()

% This function plots the SNRout values in dB for different values of the
% filter order (M), that has been precalculated and stored in a file. 
% -------------------------------------------------------------------------
% Code written by: KAUSHIK IYER 
% =========================================================================

%% ADDPATH AND LOAD FILE
addpath('../codes/Results/');
load SNR_OUT.mat;

figure();
plot(M_all,LMS_SNRout,'-d','LineWidth',4);
hold on
plot(M_all,NLMS_SNRout,'-o','LineWidth',4);
plot(M_all,ENLMS_SNRout,'-*','LineWidth',4);
plot(M_all,RLS_SNRout,'-s','LineWidth',4);
plot(M_all,KLMS_SNRout,'-+','LineWidth',4);
plot(M_all,RENLMS_SNRout,'-^','LineWidth',4);

legend('NLMS','IPNLMS','ENLMS','Hinf','KLMS','RENLMS');

%% SET TICKS
xticks(M_all);
xtickangle(0);

%% LABELS FOR THE PLOT
title('Filter Order v/s $SNR_{out}$','interpreter','latex',Fontsize=24);
xlabel('Filter Order ($M$)', 'interpreter','latex',Fontsize=24);
ylabel('$SNR_{out}$ (dB)','interpreter','latex',Fontsize=24);

%% AXIS SETTINGS 
ax = gca;
ax.FontSize = 20;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
box on
ax.LineWidth = 2;


end