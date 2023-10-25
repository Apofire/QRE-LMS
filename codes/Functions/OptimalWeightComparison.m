function OptimalWeightComparison(CleanData,NoisyData,NoisyDataIdx,LMS_settings,NLMS_settings,...
    ENLMS_settings,KLMS_settings,RLS_settings,RE_NLMS_settings,ITER)

EEGidx = NoisyDataIdx.EEG;       % Collect indices used to generate noisy dataset
d = CleanData(EEGidx(1),:)'; % Extract clean EEG signal
x = NoisyData(:,1);              % Extract noisy EEG signal

M = RLS_settings.M;  % Filter Order
Mmax = length(d);    % Max possible Filter Order

LMS_settings_max     = SetSettings(Mmax,"LMS"); 
NLMS_settings_max    = SetSettings(Mmax,"NLMS");
ENLMS_settings_max   = SetSettings(Mmax,"ENLMS");
KLMS_settings_max    = SetSettings(Mmax,"KLMS");
RLS_settings_max     = SetSettings(Mmax,"RLS");
RE_NLMS_settings_max = SetSettings(Mmax,"RE_NLMS");


% Find optimal weights using largest possible Filter order
[~,~,WoptLMS] = LMS(x,d,LMS_settings_max.w,0,ITER);

[y,~,WoptNLMS] = NLMS(x,d,NLMS_settings_max.w,NLMS_settings.mu,NLMS_settings.a,ITER);

[~,~,WoptENLMS] = ENLMS(x,d,ENLMS_settings_max.w,ENLMS_settings.mu,ENLMS_settings.a,ITER);

[~,~,WoptKLMS] = KLMS(x,d,KLMS_settings_max.lambda,KLMS_settings_max.P,...
    KLMS_settings_max.Q,KLMS_settings_max.qv,KLMS_settings_max.w,ITER);

[~,~,WoptRLS] = RLS(x,d,Mmax,ITER,RLS_settings_max.lambda,RLS_settings_max.Rinv);

[y,~,WoptRENLMS] =  RE_NLMS(x,d,RE_NLMS_settings_max.w,RE_NLMS_settings.mu,...
    RE_NLMS_settings.muNLMS,RE_NLMS_settings.a,...
    ITER,ITER);

% Run All other algorithms to get the weights
[~,~,wLMS] = LMS(x,d,LMS_settings.w,0,ITER);

[~,~,wNLMS] = NLMS(x,d,NLMS_settings.w,NLMS_settings.mu,NLMS_settings.a,ITER);

[~,~,wENLMS] = ENLMS(x,d,ENLMS_settings.w,ENLMS_settings.mu,ENLMS_settings.a,ITER);

[~,~,wKLMS] = KLMS(x,d,KLMS_settings.lambda,KLMS_settings.P,...
    KLMS_settings.Q,KLMS_settings.qv,KLMS_settings.w,ITER);

[~,~,wRLS] = RLS(x,d,RLS_settings.M, ITER, RLS_settings.lambda,RLS_settings.Rinv);

[~,~,wRENLMS] =  RE_NLMS(x,d,RE_NLMS_settings.w,RE_NLMS_settings.mu,...
    RE_NLMS_settings.muNLMS,RE_NLMS_settings.a,...
    ITER,ITER);


%% COMPUTE NORM
for j = 1:ITER
    normLMS(j)    = norm(WoptLMS(1:M,end) - wLMS(:,j));
    normNLMS(j)   = norm(WoptNLMS(1:M,end) - wNLMS(:,j));
    normENLMS(j)  = norm(WoptENLMS(1:M,end) - wENLMS(:,j));
    normKLMS(j)   = norm(WoptKLMS(1:M,end) - wKLMS(:,j));
    normRLS(j)    = norm(WoptRLS(1:M,end) - wRLS(:,j));
    normRENLMS(j) = norm(WoptRENLMS(1:M,end) - WoptRENLMS(1:M,j));
end
plot(normLMS(6:end),'LineWidth',3);
hold on
plot(flip(normNLMS(6:end)),'LineWidth',3); plot(normENLMS(6:end),'LineWidth',3);
plot(normKLMS(6:end),'LineWidth',3); plot(normRLS(6:end),'LineWidth',3);
plot(normRENLMS(6:end),'LineWidth',3);
legend('NLMS','IPNLMS','ENLMS','Hinf','KLMS','RENLMS')

title('$||\mathbf{w}_{opt} - \mathbf{w}(k)||_{2}$','interpreter','latex',Fontsize=24);
xlabel('Iteration Index (k)', 'interpreter','latex',Fontsize=24);
ylabel('Magnitude','interpreter','latex',Fontsize=24);

ax = gca;
ax.FontSize = 20;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';
box on
ax.LineWidth = 2;


end

